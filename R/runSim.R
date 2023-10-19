#' Simulations - Fit Local Models in Parallel Using GD Prior
#'
#' @param n number of subjects
#' @param L number of locations
#' @param tau number of time points
#' @param n_trn size of training set
#' @param N value chosen sufficiently larger than n to sample X from
#' @param tau_start index in vector of time points representing the start of the current subgroup
#' @param tau_end index in vector of time points representing the end of the current subgroup
#' @param modRstan model code declaring the stan model
#' @param actL number of active locations - default \eqn{actL = ceiling(0.2*L)}
#' @param actL_p_t prob of coefficient at time point t being significant for an active location - default \eqn{actL_p_{lt} = 0.8}
#' @param inactL_p_t prob of coefficient at time point t being significant for an inactive location - default \eqn{inactL_p_{lt} = 0.2}
#' @param c multiplying constant to increase the magnitude of nonzero entries in the true beta matrix - default \eqn{c = 3}
#' @param rho used to generate covariance matrix when generating data - controls degree of multicollinearity
#' @param tau0 hyperparameter used in GD prior
#' @param chains number of MCMC chains
#' @param warmup number of warmup iterations
#' @param iter total number of iterations
#' @param stan_seed seed for running stan model
#'
#' @return list containing important measures to report from simulation
#' @export
#'
#' @importFrom MASS mvrnorm
#' @importFrom rstan stan
#' @importFrom parallel detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom stats glm runif rbinom qchisq deviance predict kmeans density rnorm
#' @importFrom glmnet glmnet cv.glmnet
#' @importFrom ncvreg ncvreg cv.ncvreg
#' @importFrom foreach foreach %dopar%
#' @importFrom RADIOHEAD localFDR
#' @importFrom pROC roc coords

runSim <- function(n = 130,
                   L,
                   tau,
                   n_trn,
                   N = 2*n,
                   tau_start,
                   tau_end,
                   modRstan,
                   actL = ceiling(0.2*L),
                   actL_p_t = 0.8,
                   inactL_p_t = 0.2,
                   c = 5,
                   rho = 0.5,
                   tau0 = 10^(-5),
                   chains = 1,
                   warmup = 3000,
                   iter = 10000,
                   stan_seed = 2411
) {


    # true beta - fixed coefficient vectors and fixed set of active indices
    temp = true_beta(L = L, tau = tau, actL = ceiling(0.2*L), actL_p_t = 0.8, inactL_p_t = 0.2, c = 1)
    tempBeta = temp$beta
    beta = tempBeta * matrix(sample(c(-1, 1),
                                    size = L*tau,
                                    replace = TRUE,
                                    prob = c(0.3, 0.7)),
                             nrow = L,
                             ncol = tau)
    # beta = tempBeta
    actInd = temp$act_ind

    # binary response vector
    y = rbinom(n = n, size = 1, prob = 77/122)

    trn_ind = sample.int(n, size = n_trn) # indices of training set

    y_trn = y[trn_ind] # training set for the response vector
    y_test = y[-trn_ind] # testing set for the response vector

    # three dimensional data array for current replication
    X = genX(y = y, beta = beta, n = n, N = 2*n, L = L, tau = tau, mu = rep(0, L),
             rho = rho, S = rho^(abs(matrix(1:L, L, L) - t(matrix(1:L, L, L)))))

    # adding noise to the electrical signals (X) to account for low signal-to-noise in EEG data
    noise = rnorm(n*L*tau, mean = 0, sd = 0.25)
    X = X + array(noise, dim = dim(X))

    X_trn = X[trn_ind, , ] # training set for the data
    X_test = X[-trn_ind, , ] # testing set for the data

    modFit = localMods(y = y_trn, X = X_trn, tau_st = tau_start, tau_end = tau_end,
                       tau0 = tau0, modRstan = modRstan, chains = chains,
                       warmup = warmup, iter = iter, stan_seed = stan_seed)

    list(beta = beta,
         actInd = actInd,
         y_test = y_test,
         X_test = X_test,
         b_samps = modFit$b_samps)

    }




