\name{evolvabilityBeta}
\alias{evolvabilityBeta}

\title{
Calculate evolvability parameters along a set of selection gradients
}

\description{
\code{evolvabilityBeta} calculates (unconditional) evolvability (e), respondability (r), conditional evolvability (c), autonomy (a) and integration (i) along selection gradients given an additive-genetic variance matrix as described in Hansen and Houle (2008). 
}

\usage{
evolvabilityBeta(G, Beta, means = 1)
}

\arguments{
  \item{G}{
a variance matrix.
}
  \item{Beta}{
either a vector or a matrix of unit length selection gradients stacked column wise.
}
  \item{means}{
optional vector of trait means.
}
}

\value{
\item{Beta}{matrix of selection gradients}
\item{e}{evolvability of each selection gradient}
\item{r}{respondability of each selection gradient}
\item{c}{conditional evolvability of each selection gradient}
\item{a}{autonomy of each selection gradient}
\item{i}{integration of each selection gradient}
}

\references{
Hansen, T. F. & Houle, D. (2008) Measuring and comparing evolvability and constraint in multivariate characters. J. Evol. Biol. 21:1201-1219.
}

\author{
Geir H. Bolstad \email{geir.h.bolstad@ntnu.no}
}

\examples{
G = matrix(c(1, 1, 0, 1, 2, 1, 0, 1, 2), ncol = 3)/10
Beta = randomBeta(5, 3)
X = evolvabilityBeta(G, Beta)
summary(X)
}
\keyword{array}
\keyword{algebra}
