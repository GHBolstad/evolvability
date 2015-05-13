\name{summary.evolvabilityBeta}
\alias{summary.evolvabilityBeta}
\alias{print.summary.evolvabilityBeta}

\title{
Summarizing evolvability parameters over a set of selection gradients
}

\description{
\code{summary} method for class \code{"evolvabilityBeta"}. The returned object is suitable for printing with the \code{print.summary.evolvabilityBeta} method.
}

\usage{
\method{summary}{evolvabilityBeta}(object, ...)
}

\arguments{
  \item{object}{an object of class \code{"evolvabilityBeta"}}
  \item{...}{additional arguments affecting the summary produced.}
}

\value{
\item{Averages }{the averages of the evolvability parameters over all selection gradients}
\item{Minimum }{the minimum of the evolvability parameters over all selection gradients}
\item{Maximum }{the maximum of the evolvability parameters over all selection gradients}
}

\author{Geir H. Bolstad \email{geir.h.bolstad@ntnu.no}}

\seealso{\code{\link{evolvabilityBeta}}}

\keyword{array}
\keyword{algebra}
