\name{evolvabilityBetaMCMC}
\alias{evolvabilityBetaMCMC}

\title{
Calculate posterior distribution of evolvability parameters from a set of selection gradients
}

\description{
\code{evolvabilityBetaMCMC} calculates (unconditional) evolvability (e), respondability (r), conditional evolvability (c), autonomy (a) and integration (i) from selection gradients given the posterior distribution of an additive-genetic variance matrix. These measures and their meaning are described in Hansen and Houle (2008). 
}

\usage{
evolvabilityBetaMCMC(G_mcmc, Beta, post.dist = FALSE)
}

\arguments{
  \item{G_mcmc}{
posterior distribution of a variance matrix in the form of a table. Each row in the table must be one iteration of the posterior distribution (or bootstrap distribution). Each iteration of the matrix must be on the form as given by \code{c(x)}, where \code{x} is a matrix. A posterior distribution of a matrix in the slot \code{VCV} of a object of class \code{MCMCglmm} is by default on this form.
}
  \item{Beta}{
either a vector or a matrix of unit length selection gradients stacked column wise.
}
  \item{post.dist}{
logical: should the posterior distribution of the evolvability parameters be saved.
}
}

\value{
\item{eB }{posterior median and highest posterior density interval of evolvability for each selection gradient}
\item{rB }{posterior median and highest posterior density interval of respondability for each selection gradient}
\item{cB }{posterior median and highest posterior density interval of conditional evolvability for each selection gradient}
\item{aB }{posterior median and highest posterior density interval of autonomy for each selection gradient}
\item{iB }{posterior median and highest posterior density interval of integration for each selection gradient}
\item{Beta }{matrix of selection gradients}
\item{summary }{means of evolvability parameters across all selection gradients}
\item{post.dist }{full posterior distributions}

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

# Generating selection gradients in five random directions:
Beta = randomBeta(5, 3) 

# Calculating evolvability parameters:
x = evolvabilityBetaMCMC(G_mcmc, Beta, post.dist=TRUE)
summary(x)
}

\keyword{array}
\keyword{algebra}
\keyword{multivariate}
