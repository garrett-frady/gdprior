#' Function to Generate Data Array
#'
#' @param y binary response vector
#' @param beta coefficient matrix
#' @param n number of subjects
#' @param N value chosen sufficiently larger than n to sample X from
#' @param L number of locations
#' @param tau number of time points
#' @param mu mean of MVN distribution to generate U
#' @param rho like a correlation coefficient used to generate S
#' @param S standard deviation of MVN distribution to generate U
#'
#' @return three-dimensional array of data
#' @export

genX <- function(y,
                 beta,
                 n,
                 N,
                 L,
                 tau,
                 mu,
                 rho,
                 S
                 ) {

  # U array; entries will be generated within the simulation at each replication
  U <- array(dim = c(N, L, tau))
  # X array, i.e., the data; entries will be generated within the simulation using the U array at each replication
  X <- array(dim = c(n, L, tau))

  # generate entries for the U matrix from MVN distribution
  for (t in 1:tau) {
    for(i in 1:N) {
      U[i, , t] <- mvrnorm(1, mu, S)
    }
  }

  # compute U%*%beta and the resulting y vector. Note, the response for each...
  # ... subject across time could be different.
  # (eventually we may want to generate from a matrix distribution)
  Ubeta <- matrix(nrow = N, ncol = tau)
  for (t in 1:tau) {
    Ubeta[, t] <- crossprod(t(U[, , t]), beta[, t])
  }
  yUbeta <- matrix(data = 1, nrow = N, ncol = tau)
  yUbeta[Ubeta <= 0] <- 0

  # extract the rows of U to generate the rows of X
  for(t in 1:tau){
    ind0 <- which(yUbeta[, t] == 0)
    ind1 <- which(yUbeta[, t] == 1)
    X[y == 0, , t] <- U[sample(ind0, size = n - sum(y)), , t]
    X[y == 1, , t] <- U[sample(ind1, size = sum(y)), , t]
  }

  return(X)

}
