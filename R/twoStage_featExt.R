#' Function to perform two-stage feature extraction
#'
#' FDR in stage 1 and k-means clustering on the Riemman sums of location vectors in stage 2
#'
#' First stage of the two-stage feature extraction process
#'
#' @param b MCMC samples - dimension \eqn{Nmcmc \times L*\tau}
#' @param L number of locations
#' @param tau number of time points
#' @param stage1_thres threshold used in localFDR - defaul \eqn{0.001}
#' @param FDRc FDR level at which the FDR is controlled - default \eqn{0.05}
#'
#' @return \eqn{L \times \tau} matrix of final coefficient estimates and selected locations
#' @export

################################
# Two-Stage Feature Extraction #
################################

twoStage_featExt = function(b,
                            L,
                            tau,
                            stage1_thres = 0.001,
                            FDRc = 0.05
) {

  ######################################
  ### First-stage feature extraction ###
  ######################################

  # FDR Feature Extraction
  fdr_featExt = function(b, # MCMC samples
                         stage1_thres = stage1_thres, # threshold used in localFDR
                         FDRc = FDRc # FDR control
  ) {
    Ltau = ncol(b)
    post_inc_prob = 1 - localFDR(chain = b, thres = stage1_thres, FDRc = FDRc)
    if(sum(!is.na(post_inc_prob))>0){
      inds = which(!is.na(post_inc_prob))
      est_coef = apply(b, 2,
                       function(x){
                         obj = density(x)
                         m = obj$x[which.max(obj$y)]
                         if(m == max(obj$x[obj$x<0]) | m == min(obj$x[obj$x>0])) return(0)
                         else return(m)
                       })

      est_coef[-inds] = 0
      return(est_coef)

    } else return(rep(0, Ltau))
  }

  fdr_ests = fdr_featExt(b = b,
                         stage1_thres = stage1_thres,
                         FDRc = FDRc)

  n0b = matrix(fdr_ests, nrow = L, ncol = tau)

  #######################################
  ### Second stage feature extraction ###
  #######################################

  aucLocs = apply(X = abs(n0b), MARGIN = 1, FUN = riemSum)

  cl = kmeans(x = aucLocs, centers = 2, iter.max = 10, nstart = 25)

  actLabel = which(cl$centers == max(cl$centers))

  # indices selected by second stage of feature extraction
  bInd = which(cl$cluster == actLabel)

  list(n0b_ests = n0b,
       bInd = bInd)
}
