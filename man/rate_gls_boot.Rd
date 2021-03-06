% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rate_gls.R
\name{rate_gls_boot}
\alias{rate_gls_boot}
\title{Bootstrap of the \code{\link{rate_gls}} model fit}
\usage{
rate_gls_boot(
  object,
  n = 10,
  useLFO = TRUE,
  silent = FALSE,
  maxiter = 100,
  tol = 0.001
)
}
\arguments{
\item{object}{The output from \code{\link{rate_gls}}.}

\item{n}{The number of bootstrap samples}

\item{useLFO}{logical: when calculating the mean vector of the traits in the
'recent_evol' analysis, should the focal species be left out when
calculating the corresponding species' mean. The correct way is to use
TRUE, but in practice it has little effect and FALSE will speed up the
model fit (particularly useful when bootstrapping).}

\item{silent}{logical: whether or not the bootstrap iterations should be
printed.}

\item{maxiter}{The maximum number of iterations for updating the GLS.}

\item{tol}{tolerance for convergence. If the change in 'a' and 'b' is below
this limit between the two last iteration, convergence is reached. The
change is measured in proportion to the standard deviation of the response
for 'a' and the ratio of the standard deviation of the response to the
standard deviation of the predictor for 'b'.}
}
\value{
A list where the first slot is a table with the original estimates
  and SE from the GLS fit in the two first columns followed by the bootstrap
  estimate of the SE and the 2.5\%, 50\% and 97.5\% quantiles of the
  bootstrap distribution. The second slot contains the complete distribution.
}
\description{
\code{rate_gls_boot} performs  a parametric bootstrap of a
\code{\link{rate_gls}} model fit.
}
\examples{
# See the vignette 'Analyzing rates of evolution' and in the help 
# page of rate_gls.
}
\author{
Geir H. Bolstad
}
