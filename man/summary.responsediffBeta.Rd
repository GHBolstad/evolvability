% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/responsdiffBeta.R
\name{summary.responsediffBeta}
\alias{summary.responsediffBeta}
\title{Summarizing response differences over a set of selection gradients}
\usage{
\method{summary}{responsediffBeta}(object, ...)
}
\arguments{
\item{object}{An object of class \code{'responsediffBeta'}.}

\item{...}{Additional arguments.}
}
\value{
A list with the following components:
\tabular{llllll}{
\code{Averages} \tab\tab\tab\tab The averages of the response differences 
over all selection gradients. \cr
\code{Minimum} \tab\tab\tab\tab The minimum of the response differences 
over all selection gradients. \cr
\code{Maximum} \tab\tab\tab\tab The maximum of the response differences 
over all selection gradients.
}
}
\description{
\code{summary} method for class \code{'responsediffBeta'}.
}
\seealso{
\code{\link{responsediffBeta}}
}
\author{
Geir H. Bolstad
}
\keyword{algebra}
\keyword{array}
