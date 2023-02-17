#' Function to Run GD Simulations
#'
#' @param y response vector
#' @param X data
#' @param beta true coefficient vector
#' @param trn_ind training indices
#' @param tau0 hyperparameter used in GD prior
#' @param chains number of MCMC chains
#' @param warmup number of warmup iterations
#' @param iter total number of iterations
#' @param init_thres inital threshold in two stage feature extraction
#' @param bf_thres BF threshold in two stage feature extraction
#'
#' @return list containing important measures to report from simulation
#' @export
#'
#' @importFrom stats runif
#' @importFrom MASS mvrnorm

runSim <- function(y,
                   X,
                   beta,
                   trn_ind,
                   tau0 = 10^(-5),
                   chains = 4,
                   warmup = 3000,
                   iter = 10000,
                   init_thres = 0.0001,
                   bf_thres = 0.1

) {

  n = dim(X)[1] # number of subjects
  L = dim(X)[2] # number of locations
  tau = dim(X)[3] # number of time points

  y_trn = y[trn_ind] # training set for the response vector
  y_test = y[-trn_ind] # testing set for the response vector

  X_trn <- X[trn_ind, , ] # training set for the data
  X_test <- X[-trn_ind, , ] # testing set for the data

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

  # initial values for latent vector d and tuning parameter lambda
  d_init <- as.list(runif(n = L, min = 0, max = 10^(-5)))
  lambda_init <- as.list(runif(n = 1, min = 0, max = 10^(-3)))

  modFit = localMods(y = y_trn, X = X_trn, d_init = d_init, lambda_init = lambda_init,
                      tau0 = tau0, modRstan = modRstan, chains = chains,
                      warmup = warmup, iter = iter)

  fExt = featExt(b = modFit$b_ests, init_thres = init_thres, bf_thres = bf_thres)

  final_ests = fExt$n0b_ests
  final_ests[-c(fExt$bInd), ] = 0

  est_rates = rates_est(beta, final_ests)

  y_pred = pred(b_ests = fExt$n0b_ests,
                bInd = fExt$bInd,
                X_test = X_test)

  pred_rates = rates_pred(y_test, y_pred)

  list(modFit = modFit,
       fExt = fExt,
       est_rates = est_rates,
       pred_rates = pred_rates)

}
