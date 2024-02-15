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
                     stan_seed){

  n = dim(X)[1]
  L = dim(X)[2]
  tau = dim(X)[3]

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

                       # fit stan model
                       fit = rstan::stan(model_code = modRstan, data = dat, warmup = warmup,
                                         iter = iter, seed = stan_seed, chains = chains,
                                         cores = chains, init = init_list)

                       # mcmc estimates of beta coefficients - posterior mean over all chains
                       beta_samps = rstan::extract(fit)$beta
                       mcmc_ests = rstan::get_posterior_mean(fit, par = "beta")[, chains]

                       list(b_ests = mcmc_ests,
                            b_samps = beta_samps)

                     }

  parallel::stopCluster(cl) # Stop cluster

  # L x tau matrix of posterior mean estimates of beta - cbind vectors ...
  # ... from each time point
  b_ests = do.call(cbind, mcmc_sim[, "b_ests"])
  b_samps = do.call(cbind, mcmc_sim[, "b_samps"])


  # return the estimates of parameters obtained from the MCMC samples and ...
  # ... the MCMC samples, after building all the local models.
  list(b_ests = b_ests,
       b_samps = b_samps)
}

