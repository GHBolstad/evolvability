\name{summary.evolvabilityBetaMCMC}
\alias{summary.evolvabilityBetaMCMC}
\alias{print.summary.evolvabilityBetaMCMC}
\alias{print.evolvabilityBetaMCMC}

\title{
Summarizing posterior distribution of evolvability parameters over a set of selection gradients
}

\description{
\code{summary} method for class \code{"evolvabilityBetaMCMC"}. The returned object is suitable for printing with the \code{print.summary.evolvabilityBetaMCMC} method.
}

\usage{
\method{summary}{evolvabilityBetaMCMC}(object, ...)
}

\arguments{
  \item{object}{an object of class \code{"evolvabilityBetaMCMC"}}
  \item{...}{additional arguments affecting the summary produced.}
}

\value{
\item{Averages }{the averages of the evolvability parameters over all selection gradients}
\item{Minimum }{the minimum (given by the posterior median) of the evolvability parameters over all selection gradients}
\item{Maximum }{the maximum (given by the posterior median) of the evolvability parameters over all selection gradients}
}

\author{Geir H. Bolstad \email{geir.h.bolstad@ntnu.no}}

\seealso{\code{\link{evolvabilityBetaMCMC}}}

\keyword{array}
\keyword{algebra}
