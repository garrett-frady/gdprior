#' Generate true beta coefficients
#'
#' Active locations should display a longer streaky behavior of 1's over time
#'
#' @param L number of locations
#' @param tau number of time points
#' @param actL number of active locations - default \eqn{actL = ceiling(0.2*L)}
#' @param actL_p_t prob of coefficient at time point t being significant for an active location - default \eqn{actL_p_{lt} = 0.9}
#' @param inactL_p_t prob of coefficient at time point t being significant for an inactive location - default \eqn{inactL_p_{lt} = 0.1}
#' @param c multiplying constant to increase the magnitude of nonzero entries in the true beta matrix - default \eqn{c = 3}
#'
#' @return \eqn{L \times \tau} matrix of beta coefficients
#' @export

true_beta = function(L,
                     tau,
                     actL = ceiling(0.2*L),
                     actL_p_t = 0.9,
                     inactL_p_t = 0.1,
                     c = 3) {

  # matrix to store randomly generated beta coefficients in - L x tau
  beta = matrix(nrow = L, ncol = tau)

  # randomly select indices of active locations
  act_ind = sample.int(L, actL)

  # loop to generate time series for each location
  for (l in 1:L) {

    # if statement to determine whether the current location is active or not
    if (l %in% act_ind) {
      beta[l, ] = c*rbinom(n = tau, size = 1, prob = actL_p_t)
    } else {
      beta[l, ] = c*rbinom(n = tau, size = 1, prob = inactL_p_t)
    }
  }

  list(act_ind = act_ind,
       beta = beta)

}
