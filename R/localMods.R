#' Constructing Local Models in Parallel
#'
#' @param y binary response correspoding to \eqn{n} subjects
#' @param X \eqn{n} x \eqn{L} x \eqn{tau} - tensor (EEG data) for \eqn{n} subjects, \eqn{L} locations and \eqn{tau} time points
#' @param tau0 predetermined constant \eqn{\tau0} - default \eqn{\tau0 = 10^(-5)}
#' @param chains number of chains in MCMC algorithm - default \eqn{chains = 4}
#' @param warmup number of iterations to throw away as burnin at the start of the mcmc run - default \eqn{warmup = 3000}
#' @param iter number of iterations to run in each mcmc chain - default \eqn{iter = 10000}
#'
#' @return list containing posterior mean estimates of the coefficients and MCMC samples of the coefficients
#' @export

localMods = function(y,
                     X,
                     tau0 = 10^(-5),
                     chains = 4,
                     warmup = 3000,
                     iter = 10000){

  n = dim(X)[1]
  L = dim(X)[2]
  tau = dim(X)[3]

  # rstan MCMC code for the GD model
  modRstan = "

    functions {

      vector colSums(matrix M) {
        int ncol;
        vector[cols(M)] sums;

        ncol = cols(M);
        for (i in 1:ncol) {
          sums[i] = sum(M[, i]);
        }
        return(sums);
      }

    }

    data {

      int<lower = 0> n;
      int<lower = 0> L;
      matrix[n, L] X;
      int<lower = 0, upper = 1> y[n];
      real<lower = 10^(-10)> tau0;

    }

    transformed data {

      vector[L] tXX_diag;
      tXX_diag = colSums(X .* X);

    }

    parameters {

      vector[L] beta;
      vector<lower = 0>[L] d;
      real<lower = 0> lambda;

    }

    model {

      lambda ~ gamma(0.1, 0.2);

      for (j in 1:L) {
        d[j] ~ gamma(lambda + 0.5, tau0^2/(2*tXX_diag[j]));
      }

      for (j in 1:L) {
        beta[j] ~ normal(0, sqrt(1/d[j]));
      }

      y ~ bernoulli_logit_glm(X, 0, beta);

    }

  "

  ncores = parallel::detectCores() # detect the number of cores
  doParallel::registerDoParallel(cores = ceiling(ncores/chains)) # initiate the cluster
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
                       d_init = as.list(runif(n = L, min = 0, max = 10^(-5)))
                       lambda_init = as.list(runif(n = 1, min = 0, max = 10^(-3)))

                       # initial values for parameters; lists required by rstan
                       beta_init = as.list(glm(formula = y ~ X[, , t] + 0,
                                               family = "binomial")$coef)

                       # list of initial values; size equivalent to num of chains
                       init_list = list(one = list(beta = beta_init,
                                                   d = d_init,
                                                   lambda = lambda_init),
                                        two = list(beta = beta_init,
                                                   d = d_init,
                                                   lambda = lambda_init),
                                        three = list(beta = beta_init,
                                                     d = d_init,
                                                     lambda = lambda_init),
                                        four = list(beta = beta_init,
                                                    d = d_init,
                                                    lambda = lambda_init)
                       )

                       # fit stan model
                       fit = rstan::stan(model_code = modRstan, data = dat, warmup = warmup,
                                         iter = iter, seed = 2411, chains = chains, cores = chains,
                                         init = init_list)

                       # mcmc estimates of beta coefficients - posterior mean over all chains
                       mcmc_ests = rstan::get_posterior_mean(fit, par = "beta")[, chains + 1]
                       beta_samps = rstan::extract(fit)$beta

                       list(b_ests = mcmc_ests,
                            b_samps = beta_samps)

                     }

  # L x tau matrix of posterior mean estimates of beta - cbind vectors from each time point
  b_ests = do.call(cbind, mcmc_sim[, "b_ests"])
  b_samps = do.call(cbind, mcmc_sim[, "b_samps"])

  # return the estimates of parameters obtained from the MCMC samples...
  # ... after building all the local models.
  list(b_ests = b_ests,
       b_samps = b_samps)
}
