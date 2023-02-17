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
#' @importFrom MASS mvrnorm
#' @importFrom rstan stan get_posterior_mean extract
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom stats glm runif
#' @import foreach

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

  modFit = localMods(y = y_trn, X = X_trn, tau0 = tau0, chains = chains,
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
