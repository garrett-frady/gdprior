#' Performing two-stage feature extraction
#'
#' @param b \eqn{L} x \eqn{tau} matrix of estimates of beta
#' @param init_thres threshold in first stage to induce sparsity
#' @param bf_thres threshold to be used to compare Bayes factor
#' @param L number of locations
#' @param tau number of time points
#'
#' @return list containing nonzero estimates after first stage, active indices according to the second stage, and vector of maxBF values
#' @export

featExt = function(b,
                   init_thres,
                   bf_thres,
                   L,
                   tau) {

  ######################################
  ### First-stage feature extraction ###
  ######################################
  n0b = b
  n0b[abs(n0b) < init_thres] = 0

  #######################################
  ### Second stage feature extraction ###
  #######################################
  b_bin = matrix(1, nrow = L, ncol = tau)
  b_bin[n0b == 0] = 0

  maxBF = numeric()
  for(l in 1:L) {
    if(sum(b_bin[l, ]) == 0) {
      maxBF[l] = -10
    } else if(sum(b_bin[l, ]) == tau) {
      maxBF[l] = bf_thres
    } else maxBF[l] = bf_fn(b_bin[l, ])$max.log.BF
  }

  bInd = which(maxBF >= bf_thres)

  list(n0b_ests = n0b,
       bInd = bInd,
       maxBF = maxBF)
}
