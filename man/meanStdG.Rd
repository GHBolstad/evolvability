\name{meanStdG}
\alias{meanStdG}

\title{
Mean standardize a G-matrix
}

\description{
Mean standardization of a G-matrix.
}

\usage{
meanStdG(G, means)
}

\arguments{
  \item{G}{
a variance matrix.
}
  \item{means}{
a vector of trait means.
}
}

\value{
\code{meanStdG} returns a mean standardized variance matrix. 
}

\author{
Geir H. Bolstad
}

\examples{
G = matrix(c(1, 1, 0, 1, 4, 1, 0, 1, 2), ncol = 3)
means = c(1, 1.4, 2.1)

meanStdG(G, means)
}

\keyword{array}
\keyword{algebra}
