#' Function to Compute Error Rates of Prediction
#'
#' @param x vector containing true values of the coefficients
#' @param x_hat vector containing predicted values of the coefficients
#'
#' @return vector containing the tpr, fpr, and pe of predictions
#' @export

rates_pred <- function(x, x_hat){
  tpr = sum((x+x_hat) == 2)/sum(x)
  fpr = sum((x-x_hat) == -1)/(length(x) - sum(x))
  pe = sum((x+x_hat) == 1)/length(x)

  return(c(tpr, fpr, pe))
}
