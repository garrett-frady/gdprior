#' Function to Compute Error Rates of Estimation
#'
#' @param x # vector containing true values of the coefficients
#' @param x_hat # vector containing predicted values of the coefficients
#'
#' @return vector containing the rmse, tpr, and fpr of estimation
#' @export

rates_est <- function(x, x_hat){
  rmse = sqrt(mean((x - x_hat)^2))
  tpr = sum((x != 0 & x_hat != 0))/sum(x != 0)
  fpr = sum((x == 0 & x_hat != 0))/sum(x == 0)

  return(c(rmse, tpr, fpr))
}
