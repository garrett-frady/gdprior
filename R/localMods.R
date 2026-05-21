#' Constructing Local Models in Parallel
#'
#' @param y binary response correspoding to \eqn{n} subjects
#' @param X \eqn{n} x \eqn{L} x \eqn{tau} - tensor (EEG data) for \eqn{n} subjects, \eqn{L} locations and \eqn{tau} time points
#' @param tau0 predetermined constant \eqn{\tau0} - default \eqn{\tau0 = 10^(-5)}
#' @param modRstan model code declaring the stan model
#' @param chains number of chains in MCMC algorithm - default \eqn{chains = 4}
#' @param warmup number of iterations to throw away as burnin at the start of the mcmc run - default \eqn{warmup = 3000}
#' @param iter number of iterations to run in each mcmc chain - default \eqn{iter = 10000}
#' @param stan_seed seed for running stan model
#' @param parallelize logical; if \code{TRUE}, local models are fitted in parallel
#'   across \eqn{\tau} time points using \code{foreach \%dopar\%}. If \code{FALSE}
#'   (default), local models are fitted sequentially in a \code{for} loop. Has no
#'   effect when \eqn{\tau = 1} and \eqn{chains = 1}, which always takes the
#'   single-fit path with no loop overhead.
#'
#' @return list containing posterior mean estimates of the coefficients and MCMC samples of the coefficients
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom rstan stan
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom stats glm runif rbinom qchisq deviance predict kmeans density
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom ncvreg ncvreg cv.ncvreg
#' @importFrom foreach foreach %dopar%
#' @importFrom RADIOHEAD localFDR

localMods = function(y,
                     X,
                     tau0,
                     modRstan,
                     chains,
                     warmup,
                     iter,
                     stan_seed,
                     parallelize = FALSE){

  n = dim(X)[1]
  L = dim(X)[2]
  tau = dim(X)[3]

  # ── Internal helper: fit one Stan local model at time-slice t ──────────────
  # Shared by both tau == 1 branches to avoid code duplication.
  # Arguments are taken from the enclosing function environment.
  # Returns list(b_ests, b_samps, run_times) matching the structure produced
  # by a single row of mcmc_sim in the tau > 1 foreach branch.
  .fit_one_local <- function(t) {

    # data list passed to Stan for the local model at time point t
    dat <- list(
      n    = as.integer(n),
      L    = as.integer(L),
      X    = X[, , t],   # n x L slice at time point t
      y    = y,
      tau0 = tau0
    )

    # initialise latent scale vector d, shrinkage parameter lambda, and
    # regression coefficients beta using a plain logistic regression fit
    d_init      <- as.list(runif(n = L, min = 0, max = 1e-1))
    lambda_init <- as.list(runif(n = 1, min = 0, max = 1e-3))
    beta_init   <- as.list(glm(formula = y ~ X[, , t] + 0,
                               family  = "binomial")$coef)

    # build init list with one entry per chain (rstan requirement)
    init_list <- list()
    for (j in seq_len(chains)) {
      init_list[[paste0("chain_", j)]] <- list(beta   = beta_init,
                                               d      = d_init,
                                               lambda = lambda_init)
    }

    start_time <- Sys.time()
    fit <- rstan::stan(
      model_code = modRstan,
      data       = dat,
      warmup     = warmup,
      iter       = iter,
      seed       = stan_seed,
      chains     = chains,
      cores      = chains,   # Stan handles multi-chain parallelism internally
      init       = init_list
    )
    end_time   <- Sys.time()
    time_taken <- round(end_time - start_time, 2)

    beta_samps <- rstan::extract(fit)$beta
    mcmc_ests  <- rstan::get_posterior_mean(fit, par = "beta")[, chains]

    list(b_ests    = mcmc_ests,
         b_samps   = beta_samps,
         run_times = time_taken)
  }

  # ── Dispatch on (tau, chains) ───────────────────────────────────────────────

  if (tau == 1 && chains == 1) {

    # ── Case 1: tau = 1, chains = 1 ─────────────────────────────────────────
    # Simplest path: a single Stan fit with a single chain.
    # No cluster is created; no foreach loop is needed; no parallel overhead
    # of any kind. Stan runs sequentially on the one available time slice.

    result    <- .fit_one_local(t = 1)
    b_ests    <- matrix(result$b_ests, nrow = L, ncol = 1)
    b_samps   <- result$b_samps
    run_times <- result$run_times

  }  else if (!parallelize) {

    # ── Case 2: sequential loop over tau time points (default) ───────────────
    # Local models are fitted one at a time in a plain for loop.
    # No cluster is created; no parallel workers are spawned.
    # .fit_one_local is reused directly for each time point t.
    # This is the default path (parallelize = FALSE).

    results_list <- vector("list", tau)   # pre-allocate list of length tau

    for (t in seq_len(tau)) {
      results_list[[t]] <- .fit_one_local(t = t)
    }

    # combine results across time points into matrices, exactly as in the
    # parallel path so that downstream code is unaffected by which path ran
    b_ests    <- do.call(cbind, lapply(results_list, `[[`, "b_ests"))
    b_samps   <- do.call(cbind, lapply(results_list, `[[`, "b_samps"))
    run_times <- do.call(cbind, lapply(results_list, `[[`, "run_times"))

  } else {
    ncores = parallel::detectCores()
    cl = parallel::makeCluster(ncores/chains) # make cluster
    doParallel::registerDoParallel(cl, cores = ceiling(ncores/chains))
    # mcmc estimation in parallel over all time points
    mcmc_sim = foreach(t = 1:tau, .combine = rbind, .inorder = TRUE,
                       .errorhandling = "pass", .packages = c('rstan')) %dopar% {

                         # data list for local stan model at time point t
                         dat = list(
                           n = as.integer(n),
                           L = as.integer(L),
                           X = X[, , t], # n x L local data matrix
                           y = y,
                           tau0 = tau0
                         )

                         # initial values for latent vector d and tuning parameter lambda
                         d_init = as.list(runif(n = L, min = 0, max = 10^(-1)))
                         lambda_init = as.list(runif(n = 1, min = 0, max = 10^(-3)))

                         # initial values for parameters; lists required by rstan
                         beta_init = as.list(glm(formula = y ~ X[, , t] + 0,
                                                 family = "binomial")$coef)

                         # list of initial values; size equivalent to num of chains
                         # list of initial values; size equivalent to num of chains
                         init_list = list()
                         for(j in 1:chains) {
                           init_list[[paste0("chain_", j)]] = list(beta = beta_init,
                                                                   d = d_init,
                                                                   lambda = lambda_init)
                         }
                         # start run time for each local model
                         start_time = Sys.time()
                         # fit stan model
                         fit = rstan::stan(model_code = modRstan, data = dat, warmup = warmup,
                                           iter = iter, seed = stan_seed, chains = chains,
                                           cores = chains, init = init_list)
                         # end run time for each local model
                         end_time = Sys.time()
                         # run time to fit each local model
                         time_taken = round(end_time - start_time, 2)

                         # mcmc estimates of beta coefficients - posterior mean over all chains
                         beta_samps = rstan::extract(fit)$beta
                         mcmc_ests = rstan::get_posterior_mean(fit, par = "beta")[, chains]

                         list(b_ests = mcmc_ests,
                              b_samps = beta_samps,
                              run_times = time_taken)

                       }

    parallel::stopCluster(cl) # Stop cluster

    # L x tau matrix of posterior mean estimates of beta - cbind vectors ...
    # ... from each time point
    b_ests = do.call(cbind, mcmc_sim[, "b_ests"])
    b_samps = do.call(cbind, mcmc_sim[, "b_samps"])
    run_times = do.call(cbind, mcmc_sim[, "run_times"])
  }

  # return the estimates of parameters obtained from the MCMC samples and ...
  # ... the MCMC samples, after building all the local models.
  list(b_ests = b_ests,
       b_samps = b_samps,
       run_times = run_times)
}


