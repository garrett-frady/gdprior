#' Get coefficient estimates, prediction probabilities, and predictions from the competing methods
#'
#' @param y_trn # training set for the binary response vector; \eqn{(n - 1) \times 1} when using LOOCV
#' @param X_trn # training set for the matricized data array; \eqn{(n - 1) \times L\tau} when using LOOCV
#' @param y_test # binary response in test set; single value when using LOOCV
#' @param X_mat_test # data matrix in test set; \eqn{1 \times L\tau} when using LOOCV
#' @param tau0 # predetermined constant - default \eqn{\tau_0 = 10^(-5)}
#'
#' @return list of coefficient estimates, prediction probabilities, and predicted responses for competing methods
#' @export

getCompMethRes = function(y_trn, # binary response vector - training set
                          X_trn, # data array - training set
                          y_test, # binary response vector - testing set
                          X_mat_test, # data array - testing set
                          tau0 = 10^(-5) # predetermined constant
){

  ##### data setup #####
  tau = dim(X_trn)[3] # number of time points
  L = dim(X_trn)[2] # number of locations
  n_trn = dim(X_trn)[1] # number of subjects

  X_mat_trn = matrix(X_trn, nrow = n_trn, ncol = L*tau)

  # Ridge Estimator - CV
  ridge_cv_model <- cv.glmnet(X_mat_trn, y_trn, intercept=FALSE, family="binomial", alpha=0)
  ridge_mod <- glmnet(X_mat_trn, y_trn, lambda = ridge_cv_model$lambda.min, intercept=FALSE, family="binomial", alpha=0)
  beta_ridge_cv <- ridge_mod$beta


  # Comparison: Lasso (CV based)
  cv_lasso_model <- cv.glmnet(X_mat_trn, y_trn, family = "binomial", alpha = 1, intercept = FALSE)
  model_lasso_cv <- glmnet(X_mat_trn, y_trn, family = "binomial", alpha = 1, intercept = FALSE, lambda = cv_lasso_model$lambda.min)
  beta_lasso_cv <- model_lasso_cv$beta
  pred_lasso_cv <- as.numeric(predict(object = model_lasso_cv,
                                      type = "class",
                                      newx = X_mat_test))
  probs_lasso_cv <- as.numeric(predict(object = model_lasso_cv,
                                       type = "response",
                                       newx = X_mat_test))


  # Comparison: Adaptive Lasso with Ridge Weights (CV based)
  hat_W <- matrix(1/abs(beta_ridge_cv), n_trn, L*tau, byrow = TRUE)
  xtilde <- X_mat_trn/hat_W
  alasso_cv_model <- cv.glmnet(xtilde, y_trn, intercept = FALSE, family = "binomial", alpha = 1)
  model_alasso_cv <- glmnet(xtilde, y_trn, lambda = alasso_cv_model$lambda.min, intercept = FALSE, family = "binomial", alpha = 1)
  beta_tilde_cv <- model_alasso_cv$beta
  # Transform beta
  beta_alasso_cv <- beta_tilde_cv/hat_W[1, ]
  pred_alasso_cv <- as.numeric(predict(object = model_alasso_cv,
                                       type = "class",
                                       newx = X_mat_test))

  probs_alasso_cv <- as.numeric(predict(object = model_alasso_cv,
                                        type = "response",
                                        newx = X_mat_test))


  # Comparison: Elastic Net (CV based)
  enet_cv_model <- cv.glmnet(X_mat_trn, y_trn, intercept = FALSE, family = "binomial", alpha = 0.5)
  model_enet_cv <- glmnet(X_mat_trn, y_trn, lambda = enet_cv_model$lambda.min, intercept = FALSE, family = "binomial", alpha = 0.5)
  beta_enet_cv <- model_enet_cv$beta
  pred_enet_cv <- as.numeric(predict(object = model_enet_cv,
                                     type = "class",
                                     newx = X_mat_test))

  probs_enet_cv <- as.numeric(predict(object = model_enet_cv,
                                      type = "response",
                                      newx = X_mat_test))


  # Comparison: MCP (CV based)
  mcp_cv_model <- cv.ncvreg(X_mat_trn, y_trn, family = "binomial", penalty = "MCP", alpha = 1, intercept = FALSE)
  min_lambda_mcp <- mcp_cv_model$lambda.min
  model_mcp_cv <- glmnet(X_mat_trn, y_trn, family = "binomial", alpha = 1, intercept = FALSE, lambda = min_lambda_mcp)
  beta_mcp_cv <- model_mcp_cv$beta
  pred_mcp_cv <- as.numeric(predict(object = model_mcp_cv,
                                    type = "class",
                                    newx = X_mat_test))

  probs_mcp_cv <- as.numeric(predict(object = model_mcp_cv,
                                     type = "response",
                                     newx = X_mat_test))


  # Comparison: SCAD (CV based)
  scad_cv_model <- cv.ncvreg(X_mat_trn, y_trn, family = "binomial", penalty = "SCAD", alpha = 1, intercept = FALSE)
  min_lambda_scad <- scad_cv_model$lambda.min
  model_scad_cv <- glmnet(X_mat_trn, y_trn, family = "binomial", alpha = 1, intercept = FALSE, lambda = min_lambda_scad)
  beta_scad_cv <- model_scad_cv$beta
  pred_scad_cv <- as.numeric(predict(object = model_scad_cv,
                                     type = "class",
                                     newx = X_mat_test))

  probs_scad_cv <- as.numeric(predict(object = model_scad_cv,
                                      type = "response",
                                      newx = X_mat_test))


  list(beta_lasso_cv = beta_lasso_cv,
       beta_alasso_cv = beta_alasso_cv,
       beta_enet_cv, beta_enet_cv,
       beta_mcp_cv = beta_mcp_cv,
       beta_scad_cv = beta_scad_cv,
       pred_lasso_cv = pred_lasso_cv,
       pred_alasso_cv = pred_alasso_cv,
       pred_enet_cv = pred_enet_cv,
       pred_mcp_cv = pred_mcp_cv,
       pred_scad_cv = pred_scad_cv,
       probs_lasso_cv = probs_lasso_cv,
       probs_alasso_cv = probs_alasso_cv,
       probs_enet_cv = probs_enet_cv,
       probs_mcp_cv = probs_mcp_cv,
       probs_scad_cv = probs_scad_cv)

}


