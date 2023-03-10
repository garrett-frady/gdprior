% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runSim.R
\name{runSim}
\alias{runSim}
\title{Function to Run GD Simulations}
\usage{
runSim(
  n,
  L,
  tau,
  n_trn,
  seed = 1234,
  ncores = parallel::detectCores(),
  N = 2 * n,
  actL = ceiling(0.2 * L),
  actL_p_t = 0.9,
  inactL_p_t = 0.1,
  c = 5,
  mu = rep(0, L),
  rho = 0.5,
  S = rho^(abs(matrix(1:L, L, L) - t(matrix(1:L, L, L)))),
  tau0 = 10^(-5),
  chains = 1,
  warmup = 3000,
  iter = 10000,
  init_thres = 10^(-4),
  bf_thres = 0.1,
  stan_seed = 2411,
  log.K = seq(0, 5, 0.1)
)
}
\arguments{
\item{n}{number of subjects}

\item{L}{number of locations}

\item{tau}{number of time points}

\item{n_trn}{size of training set}

\item{seed}{seed for reproducibility}

\item{ncores}{number of cores to use}

\item{N}{value chosen sufficiently larger than n to sample X from}

\item{actL}{number of active locations - default \eqn{actL = ceiling(0.2*L)}}

\item{actL_p_t}{prob of coefficient at time point t being significant for an active location - default \eqn{actL_p_{lt} = 0.9}}

\item{inactL_p_t}{prob of coefficient at time point t being significant for an inactive location - default \eqn{inactL_p_{lt} = 0.1}}

\item{c}{multiplying constant to increase the magnitude of nonzero entries in the true beta matrix - default \eqn{c = 3}}

\item{mu}{mean of MVN distribution to generate U}

\item{rho}{like a correlation coefficient used to generate S}

\item{S}{standard deviation of MVN distribution to generate U}

\item{tau0}{hyperparameter used in GD prior}

\item{chains}{number of MCMC chains}

\item{warmup}{number of warmup iterations}

\item{iter}{total number of iterations}

\item{init_thres}{inital threshold in two stage feature extraction}

\item{bf_thres}{BF threshold in two stage feature extraction}

\item{stan_seed}{seed for running stan model}

\item{log.K}{values to consider for K, parameter in beta distribution for the passive model}
}
\value{
list containing important measures to report from simulation
}
\description{
Function to Run GD Simulations
}
