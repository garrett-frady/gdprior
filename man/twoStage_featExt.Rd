% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twoStage_featExt.R
\name{twoStage_featExt}
\alias{twoStage_featExt}
\title{Function to perform two-stage feature extraction}
\usage{
twoStage_featExt(b, L, tau, stage1_thres = 0.001, FDRc = 0.05)
}
\arguments{
\item{b}{MCMC samples - dimension \eqn{Nmcmc \times L*\tau}}

\item{L}{number of locations}

\item{tau}{number of time points}

\item{stage1_thres}{threshold used in localFDR - defaul \eqn{0.001}}

\item{FDRc}{FDR level at which the FDR is controlled - default \eqn{0.05}}
}
\value{
\eqn{L \times \tau} matrix of final coefficient estimates and selected locations
}
\description{
FDR in stage 1 and k-means clustering on the Riemman sums of location vectors in stage 2
}
\details{
First stage of the two-stage feature extraction process
}
