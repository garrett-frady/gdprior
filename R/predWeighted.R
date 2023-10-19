#' Weighted Subject-Level Prediction
#'
#' Component-wise products of the local weights and local prediction probabilities
#'
#' @param b_ests \eqn{L \times \tau} matrix of coefficient estimates from two-stage feature extraction
#' @param bInd indices of selected locations from two-stage feature extraction
#' @param X_test \eqn{L \times \tau} dimensional matrix of test data for subject i
#' @param thres threshold used to make final prediction for subject i
#'
#' @return list containing subject-level predictions, local probabilities, local weights, and weighted predictions
#' @export

############################
# Subject-Level Prediction #
############################

predWeighted = function(b_ests, # estimates of beta returned from two stage,
                        bInd, # indices of selected locations from two-stage feature extraction
                        X_test, # test data array
                        thres # threshold used for final prediction
){

  tau = ifelse(is.null(dim(X_test)), 1, dim(X_test)[2])

  # inverse logit link function - standard logistic function
  h <- function(x) {
    return(1/(1 + exp(-x)))
  }

  if (tau < 2) {
    if (length(bInd) == 1) {
      Xb = X_test[bInd]*b_ests[bInd]
    } else {
      Xb = X_test[bInd]%*%b_ests[bInd]
    }
  } else {
    if (length(bInd) == 1) {
      Xb = sapply(1:tau, function(t) X_test[bInd, t]*b_ests[bInd, t])
    } else {
      Xb = sapply(1:tau, function(t) X_test[bInd, t]%*%b_ests[bInd, t])
    }
  }

  # local prediction probabilities over all time points; for subject i in ...
  # ... the test set (the only subject in the test set when using LOOCV)
  p_it = c()
  for (t in 1:tau) {
    p_it[t] = h(Xb[t])
  }

  # denominator used to compute the local weights
  denWeights = sum((p_it - 0.5)^2)

  # local weights subject i over all time points
  w_it = ((p_it - 0.5)^2)/denWeights

  # weighted sum over local prediction probabilities for subject i
  weightedPred = sum(w_it*p_it)

  # final predicted binary response for subject i
  y_p = ifelse(weightedPred > thres, 1, 0)

  list(y_p = y_p,
       weightedPred = weightedPred,
       p_it = p_it,
       w_it = w_it,
       Xb = Xb)
}
