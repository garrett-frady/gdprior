% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getCompMethRes.R
\name{getCompMethRes}
\alias{getCompMethRes}
\title{Get coefficient estimates, prediction probabilities, and predictions from the competing methods}
\usage{
getCompMethRes(y_trn, X_trn, y_test, X_mat_test, tau0 = 10^(-5))
}
\arguments{
\item{y_trn}{# training set for the binary response vector; \eqn{(n - 1) \times 1} when using LOOCV}

\item{X_trn}{# training set for the matricized data array; \eqn{(n - 1) \times L\tau} when using LOOCV}

\item{y_test}{# binary response in test set; single value when using LOOCV}

\item{X_mat_test}{# data matrix in test set; \eqn{1 \times L\tau} when using LOOCV}

\item{tau0}{# predetermined constant - default \eqn{\tau_0 = 10^(-5)}}
}
\value{
list of coefficient estimates, prediction probabilities, and predicted responses for competing methods
}
\description{
Get coefficient estimates, prediction probabilities, and predictions from the competing methods
}
