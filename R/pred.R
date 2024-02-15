#' Subject Level Prediction Using Length of Longest Run
#'
#' @param b_ests \eqn{L \times \tau} matrix of coefficient estimates from two-stage feature extraction
#' @param bInd indices of selected locations from two-stage feature extraction
#' @param X_test \eqn{m \times L \times \tau} dimensional array of test data
#'
#' @return binary predicted response vector of length \eqn{m}
#' @export

pred = function(b_ests, # estimates of beta returned from two stage,
                bInd, # indices of selected locations from two-stage feature extraction
                X_test # test data array
){

  m = dim(X_test)[1]
  tau = dim(X_test)[3]

  # logit link function
  h <- function(x) {
    return(1/(1 + exp(-x)))
  }

  if(length(bInd) == 1){
    Xb = sapply(1:tau, function(t) X_test[, bInd, t]*b_ests[bInd, t])
  } else{
    Xb = sapply(1:tau, function(t) X_test[, bInd, t]%*%b_ests[bInd, t])
  }

  temp = matrix(nrow = m, ncol = tau)
  for(t in 1:tau) temp[,t] = rbinom(m, 1, h(Xb[,t]))
  y_p = sapply(1:m, function(i){
    g = rle(temp[i,])
    g$values[which.max(g$lengths)]
  })

  y_p
}
