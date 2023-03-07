#' Function to Run GD Simulations
#'
#' @param n number of subjects
#' @param L number of locations
#' @param tau number of time points
#' @param n_trn size of training set
#' @param seed seed for reproducibility
#' @param ncores number of cores to use
#' @param N value chosen sufficiently larger than n to sample X from
#' @param actL number of active locations - default \eqn{actL = ceiling(0.2*L)}
#' @param actL_p_t prob of coefficient at time point t being significant for an active location - default \eqn{actL_p_{lt} = 0.9}
#' @param inactL_p_t prob of coefficient at time point t being significant for an inactive location - default \eqn{inactL_p_{lt} = 0.1}
#' @param c multiplying constant to increase the magnitude of nonzero entries in the true beta matrix - default \eqn{c = 3}
#' @param mu mean of MVN distribution to generate U
#' @param rho like a correlation coefficient used to generate S
#' @param S standard deviation of MVN distribution to generate U
#' @param tau0 hyperparameter used in GD prior
#' @param chains number of MCMC chains
#' @param warmup number of warmup iterations
#' @param iter total number of iterations
#' @param init_thres inital threshold in two stage feature extraction
#' @param bf_thres BF threshold in two stage feature extraction
#' @param stan_seed seed for running stan model
#'
#' @return list containing important measures to report from simulation
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom rstan stan get_posterior_mean extract
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom stats glm runif rbinom
#' @importFrom LearnBayes laplace
#' @import foreach

runSim <- function(n,
                   L,
                   tau,
                   n_trn,
                   seed = 1234,
                   ncores = parallel::detectCores(),
                   N = 2*n,
                   actL = ceiling(0.2*L),
                   actL_p_t = 0.9,
                   inactL_p_t = 0.1,
                   c = 3,
                   mu = rep(0, L),
                   rho = 0.5,
                   S = rho^(abs(matrix(1:L, L, L) - t(matrix(1:L, L, L)))),
                   tau0 = 10^(-5),
                   chains = 1,
                   warmup = 3000,
                   iter = 10000,
                   init_thres = 0.0001,
                   bf_thres = 0.1,
                   stan_seed = 2411

) {

  # set seed for reprodicibility
  set.seed(seed)
  # generate true coefficient matrix - L x tau
  gen_beta = true_beta(L = L, tau = tau, actL = actL, actL_p_t = actL_p_t,
                       inactL_p_t = inactL_p_t, c = c)
  beta = gen_beta$beta # true beta coefficients
  act_ind = gen_beta$act_ind # indices of true active locations

  y = rbinom(n = n, size = 1, prob = 77/122)

  X = genX(y = y, beta = beta, n = n, N = N, L = L, tau = tau, mu = mu, rho = rho, S = S)

  trn_ind = sample.int(n, size = n_trn) # indices of training set

  y_trn = y[trn_ind] # training set for the response vector
  y_test = y[-trn_ind] # testing set for the response vector

  X_trn <- X[trn_ind, , ] # training set for the data
  X_test <- X[-trn_ind, , ] # testing set for the data

  modFit = localMods(y = y_trn, X = X_trn, n_trn = n_trn, L = L, tau = tau,
                     tau0 = tau0, chains = chains, ncores = ncores,
                     warmup = warmup, iter = iter, stan_seed = stan_seed)

  fExt = featExt(b = modFit$b_ests, init_thres = init_thres, bf_thres = bf_thres,
                 L = L, tau = tau)

  final_ests = fExt$n0b_ests
  final_ests[-c(fExt$bInd), ] = 0

  est_rates = rates_est(beta, final_ests)

  y_pred = pred(b_ests = fExt$n0b_ests,
                bInd = fExt$bInd,
                X_test = X_test)

  pred_rates = rates_pred(y_test, y_pred)

  list(beta = beta,
       act_ind = act_ind,
       beta_ests = modFit$b_ests,
       beta_samps = modFit$b_samps,
       n0beta_ests = fExt$n0b_ests,
       bInd = fExt$bInd,
       est_rates = est_rates,
       pred_rates = pred_rates)

}




