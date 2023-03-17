#' Title
#'
#' @param y binary response vector corresponding to \eqn{n} subjects
#' @param X \eqn{n \times L \times \tau} dimensional tensor of EEG data for \eqn{n} subjects, \eqn{L} locations, and \eqn{\tau} time points
#' @param tau0 hyperparameter used in GD prior
#' @param n_trn size of training set
#' @param alpha significance level used to compute chi-squared threshold in second stage; controls selection rate
#' @param stage1_thres threshold in first stage to induce sparsity; fixed quantity chosen sufficiently small
#' @param ncores number of cores to use when fitting local models in parallel
#' @param chains number of chains to use in MCMC sampling
#' @param warmup number of warmup iterations to use in MCMC sampling
#' @param iter number of total iterations in MCMC sampling
#' @param seed random seed for reproducibility
#' @param stan_seed random seed used in stan code
#'
#' @return list containing important measures to report from the GD fitted model
#' @export

runEEG <- function(y,
                   X,
                   tau0 = 10^(-5),
                   n_trn = 100,
                   alpha = 0.05,
                   stage1_thres = 10^(-3),
                   ncores = detectCores(),
                   chains = 1,
                   warmup = 300,
                   iter = 1000,
                   seed = 1234,
                   stan_seed = 2411) {

  tau = dim(X)[1] # number of time points
  L = dim(X)[2] # number of locations
  n = dim(X)[3] # number of subjects

  # training and testing set
  trn_ind = sample.int(n, size = n_trn) # indices of training set
  X_trn = X[trn_ind, , ] # training set - data
  y_trn = y[trn_ind] # training set - response
  X_test = X[-trn_ind, , ] # testing set - data
  y_test = y[-trn_ind] # testing set - response

  # fit the local models and conduct MCMC estimation in parallel
  modFit = localMods(y = y_trn, X = X_trn, n_trn = n_trn, L = L, tau = tau,
                     tau0 = tau0, ncores = ncores, chains = chains,
                     warmup = warmup, iter = iter, stan_seed = stan_seed)

  # threshold for second stage - chi-squared upper (1 - alpha)-quantile
  stage2_thres = qchisq(p = 1 - alpha, df = tau)

  # two stage feature extraction
  fExt = featExt(b = modFit$b_ests, d = modFit$d_ests,
                 L = L, tau = tau, alpha = alpha,
                 stage1_thres = stage1_thres, stage2_thres = stage2_thres)

  act_b_samps = modFit$b_samps[, (1:L*tau) %% L %in% fExt$bInd]

  final_ests = fExt$n0b_ests
  final_ests[-c(fExt$bInd), ] = 0

  est_rates = rates_est(beta, final_ests)

  y_pred = pred(b_ests = fExt$n0b_ests,
                bInd = fExt$bInd,
                X_test = X_test)

  pred_rates = rates_pred(y_test, y_pred)

  list(beta_ests = modFit$b_ests,
       act_beta_samps = act_b_samps,
       bInd = fExt$bInd,
       est_rates = est_rates,
       pred_rates = pred_rates)

}
