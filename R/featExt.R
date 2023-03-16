#' Performing two-stage feature extraction
#'
#' @param b \eqn{L \times \tau} matrix of estimates of beta; posterior means of MCMC samples
#' @param d \eqn{L \times \tau} matrix of estimates of auxiliary variable d; posterior means of MCMC samples
#' @param L number of locations
#' @param tau number of time points
#' @param alpha significance level used to compute chi-squared threshold in second stage; controls selection rate
#' @param stage1_thres threshold in first stage to induce sparsity; fixed quantity chosen sufficiently small
#' @param stage2_thres threshold in second stage; chi-squared upper (1 - alpha)-quantile
#'
#' @return list containing nonzero estimates after first stage, active indices according to the second stage

featExt = function(b,
                   d,
                   L,
                   tau,
                   alpha,
                   stage1_thres,
                   stage2_thres) {

  ######################################
  ### First-stage feature extraction ###
  ######################################
  n0b = b
  n0b[abs(n0b) < stage1_thres] = 0

  #######################################
  ### Second stage feature extraction ###
  #######################################
  # ind.N(0, 1) r.v's; ind. time points over a given location - we can say ...
  # ,,, this as we built independent local models at each time point
  bTilde = sqrt(d)*b
  # ind. chi-squared r.v.'s over all locations
  bChisq = sapply(1:L, function(l) sum(bTilde[l, ]^2))

  # indices selected by second stage of feature extraction
  bInd = which(bChisq > stage2_thres)

  list(n0b_ests = n0b,
       bInd = bInd)
}
