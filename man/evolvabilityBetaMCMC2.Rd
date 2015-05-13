\name{evolvabilityBetaMCMC2}
\alias{evolvabilityBetaMCMC2}
\alias{print.evolvabilityBetaMCMC2}

\title{
Calculate posterior distribution of evolvability parameters from a selection gradient with uncertainty
}

\description{
\code{evolvabilityBetaMCMC2} calculates (unconditional) evolvability (e), respondability (r), conditional evolvability (c), autonomy (a) and integration (i) along a selection gradient with uncertainty.
}

\usage{
evolvabilityBetaMCMC2(G_mcmc, Beta_mcmc, post.dist = FALSE)
}

\arguments{
  \item{G_mcmc}{
posterior distribution of a variance matrix in the form of a table. Each row in the table must be one iteration of the posterior distribution (or bootstrap distribution). Each iteration of the matrix must be on the form as given by \code{c(x)}, where \code{x} is a matrix. A posterior distribution of a matrix in the slot \code{VCV} of a object of class \code{MCMCglmm} is by default on this form.
}
  \item{Beta_mcmc}{
posterior distribution of a unit length selection gradient where iterations are given row wise.
}
  \item{post.dist}{
logical: should the posterior distribution of the evolvability parameters be saved.
}
}

\value{
\item{Beta.median}{posterior median and highest posterior density interval of the selection gradient}
\item{summary}{posterior median and highest posterior density interval of evolvability parameters}
\item{post.dist}{full posterior distributions for the evolvability parameters}
}

\references{
Hansen, T. F. & Houle, D. (2008) Measuring and comparing evolvability and constraint in multivariate characters. J. Evol. Biol. 21:1201-1219.
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

# Simulating a posterior distribution (or bootstrap distribution) 
# of a unit length selection gradient:
Beta = randomBeta(1, 3)
Beta.mcmc = sapply(c(Beta), function(x) rnorm(10, x, 0.01))
Beta.mcmc = t(apply(Beta.mcmc, 1, function(x) x/sqrt(sum(x^2))))

# Running the model:
evolvabilityBetaMCMC2(G_mcmc, Beta_mcmc = Beta.mcmc, post.dist=TRUE)
}

\keyword{array}
\keyword{algebra}
\keyword{multivariate}
