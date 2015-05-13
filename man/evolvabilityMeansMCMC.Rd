\name{evolvabilityMeansMCMC}
\alias{evolvabilityMeansMCMC}
\alias{print.evolvabilityMeansMCMC}

\title{
Calculate posterior distribution of average evolvability parameters of a G-matrix
}

\description{
\code{evolvabilityMeans} calculates the average (unconditional) evolvability (e), respondability (r), conditional evolvability (c), autonomy (a) and integration (i) given the posterior distribution of a additive-genetic variance matrix using the approximation formulas described in Hansen and Houle (2008, 2009).
}

\usage{
evolvabilityMeansMCMC(G_mcmc)
}

\arguments{
  \item{G_mcmc}{
posterior distribution of a variance matrix in the form of a table. Each row in the table must be one iteration of the posterior distribution (or bootstrap distribution). Each iteration of the matrix must be on the form as given by \code{c(x)}, where \code{x} is a matrix. A posterior distribution of a matrix in the slot \code{VCV} of a object of class \code{MCMCglmm} is by default on this form.
}
}

\details{The formulas for calculating the evolvability parameters are approximations, except for the formula for unconditional evolvability which is exact. The bias of the approximations depends on the dimensionality of the G-matrix, with higher bias for few dimensions (see Hansen and Houle 2008). In most instances we recommend to estimate the averages of the evolvability parameters using \code{evolavbilityBetaMCMC} on a large number of random selection gradients. The maximum and minimum evolvability, which are also the maximum and minimum respondability and conditional evolvability, equals the largest and smallest eigenvalue of the G-matrix, respectively.}

\value{
\item{post.dist}{posterior distribution of the average evolvability parameters}
\item{post.medians}{posterior medians and HPD interval of the average evolvability parameters}
}


\references{
Hansen, T. F. & Houle, D. (2008) Measuring and comparing evolvability and constraint in multivariate characters. J. Evol. Biol. 21:1201-1219.
Hansen, T. F. & Houle, D. (2009) Corrigendum. J. Evol. Biol. 22:913-915.
}

\author{
Geir H. Bolstad \email{geir.h.bolstad@ntnu.no}
}

\examples{
# Simulating a posterior distribution 
# (or bootstrap distribution) of a G-matrix:
G = matrix(c(1, 1, 0, 1, 4, 1, 0, 1, 2), ncol = 3)
G_mcmc = sapply(c(G), function(x) rnorm(10, x, 0.01)) 
G_mcmc = t(apply(G_mcmc, 1, function(x){
  G = matrix(x, ncol=sqrt(length(x)))
  G[lower.tri(G)] = t(G)[lower.tri(G)]
  c(G)
}))

# Simulating a posterior distribution 
# (or bootstrap distribution) of trait means:
means = c(1, 1.4, 2.1)
means_mcmc = sapply(means, function(x) rnorm(10, x, 0.01)) 

# Mean standardizing the G-matrix:
G_mcmc = meanStdGMCMC(G_mcmc, means_mcmc)

# Estimating average evolvability paramters:
evolvabilityMeansMCMC(G_mcmc)
}

\keyword{array}
\keyword{algebra}
