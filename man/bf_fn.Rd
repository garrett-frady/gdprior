% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bf_fn.R
\name{bf_fn}
\alias{bf_fn}
\title{Bayes Factor Feature Extraction Approach}
\usage{
bf_fn(bin_vec, log.K)
}
\arguments{
\item{bin_vec}{vector of 0's and 1's}

\item{log.K}{values to consider for K, parameter in beta distribution for the passive model}
}
\value{
list of Bayes factor statistics
}
\description{
Second stage of the two-stage feature extraction process
}
\details{
Slightly adjusted version of the bayes_factor_logK function from BayesTestStreak
}
