#' Process MCMC samples for an individual replication
#'
#' Compute estimation rates, proportion of correct and incorrect indices, prediction error rates, and AUC
#'
#' @param beta # true coefficient vector
#' @param b_samps # MCMC samples obtained from GD model
#' @param y_test # response vector - test set
#' @param X_test # data array - test set
#' @param actInd # indices of predetermined active locations
#' @param stage1_thres threshold in first stage to induce sparsity; fixed quantity chosen sufficiently small
#' @param FDRc FDR level at which the FDR is controlled - default \eqn{0.05}
#'
#' @return coefficient estimates, selected locations, estimation and prediction error rates, proportion of correct and incorrect indices, subject-level predictions, ROC, AUC, best threshold
#' @export

procSim <- function(beta,
                    b_samps,
                    y_test,
                    X_test,
                    actInd,
                    stage1_thres = 10^(-3),
                    FDRc = 0.05) {

  L = dim(X_test)[2]
  tau = dim(X_test)[3]
  # feature extraction
  fExt = twoStage_featExt(b = b_samps, L = L, tau = tau, stage1_thres = stage1_thres, FDRc = FDRc)

  # store selected locations and estimates from two-stage feature extraction
  bInd = fExt$bInd
  n0b = fExt$n0b_ests

  estRates = rates_est(beta, n0b)

  # proportion of correctly selected indices and incorrectly selected indices
  corr_indices = length(which(bInd %in% actInd))/length(actInd)
  incorr_indices = length(which(!(bInd %in% actInd)))/length(actInd)

  # final predictions
  predObj = predWeighted(b_ests = n0b,
                         bInd = bInd,
                         X_test = X_test,
                         thres = 0.5)
  wPred = predObj$weightedPred
  yp = predObj$y_p
  p_it = predObj$p_it
  w_it = predObj$w_it

  # prediction error rates
  predRates = rates_pred(y_test, yp)

  # ROC and AUC approach
  ROC = pROC::roc(y_test, wPred)
  AUC = ROC$auc
  bestThres = pROC::coords(ROC, "best", ret = "threshold")

  list(n0b = n0b,
       bInd = bInd,
       estRates = estRates,
       corr_indices = corr_indices,
       incorr_indices = incorr_indices,
       yp = yp,
       predRates = predRates,
       ROC = ROC,
       AUC = AUC,
       bestThres = bestThres)

}
