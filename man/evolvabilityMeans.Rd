\name{evolvabilityMeans}
\alias{evolvabilityMeans}

\title{
Calculate average evolvability parameters of a G-matrix
}

\description{
\code{evolvabilityMeans} calculates the average (unconditional) evolvability (e), respondability (r), conditional evolvability (c), autonomy (a) and integration (i) of a additive-genetic variance matrix using the approximation formulas described in Hansen and Houle (2008, 2009). 
}

\usage{
evolvabilityMeans(G, means = 1)
}

\arguments{
  \item{G}{
a variance matrix.
}
  \item{means}{
optional vector of trait means, for mean standardization.
}
}

\details{The formulas for calculating the evolvability parameters are approximations, except for the formula for unconditional evolvability which is exact. The bias of the approximations depends on the dimensionality of the G-matrix, with higher bias for few dimensions (see Hansen and Houle 2008). In most instances we recommend to estimate the averages of the evolvability parameters using \code{evolavbilityBeta} on a large number of random selection gradients. The maximum and minimum evolvability, which are also the maximum and minimum respondability and conditional evolvability, equals the largest and smallest eigenvalue of the G-matrix, respectively.}

\value{
\item{e_mean}{average (unconditional) evolvability}
\item{e_min}{minimum evolvability}
\item{e_max}{maximum evolvability}
\item{r_mean}{average respondability}
\item{c_mean}{average conditional evolvability}
\item{a_mean}{average autonomy}
\item{i_mean}{average integration}
}


\references{
Hansen, T. F. & Houle, D. (2008) Measuring and comparing evolvability and constraint in multivariate characters. J. Evol. Biol. 21:1201-1219.
Hansen, T. F. & Houle, D. (2009) Corrigendum. J. Evol. Biol. 22:913-915.
}

\author{
Geir H. Bolstad \email{geir.h.bolstad@ntnu.no}
}

\examples{
G = matrix(c(1, 1, 0, 1, 2, 1, 0, 1, 2), ncol = 3)
evolvabilityMeans(G)
}
\keyword{array}
\keyword{algebra}
