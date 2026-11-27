#' Calculate posterior distribution of average evolvability parameters of a
#' G-matrix
#'
#' \code{evolvabilityMeansMCMC} calculates the average (unconditional) evolvability
#' (e), respondability (r), conditional evolvability (c), autonomy (a) and
#' integration (i) given the posterior distribution of a additive-genetic
#' variance matrix using the approximation formulas described in Hansen and
#' Houle (2008, 2009).
#'
#' @param G_mcmc the posterior distribution of a variance matrix in the form of
#'   a table. Each row in the table must be one iteration of the posterior
#'   distribution (or bootstrap distribution). Each iteration of the matrix must
#'   be on the form as given by \code{c(x)}, where \code{x} is a matrix. A
#'   posterior distribution of a matrix in the slot \code{VCV} of a object of
#'   class \code{MCMCglmm} is by default on this form.
#' @details The equations for calculating the evolvability parameters are
#'   approximations, except for the minimum, maximum and unconditional
#'   evolvability which are exact. The bias of the approximations depends on the
#'   dimensionality of the G-matrix, with higher bias for few dimensions (see
#'   Hansen and Houle 2008). For low dimensional G-matrices, we recommend
#'   estimating the averages of the evolvability parameters using
#'   \code{evolvabilityBetaMCMC} over many random selection gradients (
#'   \code{randomBeta}). The maximum and minimum evolvability, which
#'   are also the maximum and minimum respondability and conditional
#'   evolvability, equals the largest and smallest eigenvalue of the G-matrix,
#'   respectively. An error is given if any draw of \code{G_mcmc} is not
#'   symmetric. If \code{solve(G)} fails for a draw because the matrix is
#'   singular or computationally singular, a warning is issued, and the
#'   conditional evolvability, autonomy and integration are \code{NA} for that
#'   draw and in their posterior summaries; the other parameters are still
#'   calculated. A warning is also issued if draws are not positive definite
#'   (they have negative eigenvalues); the parameters are then still calculated
#'   but may not be meaningful for those draws.
#' @return An object of \code{class} \code{'evolvabilityMeansMCMC'}, which is a
#'   list with the following components:
#' \tabular{llllll}{
#' \code{post.dist} \tab\tab\tab\tab The posterior distribution of the average
#' evolvability parameters. \cr
#' \code{post.medians} \tab\tab\tab\tab The posterior medians and HPD interval
#' of the average evolvability parameters.
#' }
#' @references
#' Hansen, T. F. & Houle, D. (2008) Measuring and comparing evolvability and
#'  constraint in multivariate characters. J. Evol. Biol. 21:1201-1219. \cr
#' Hansen, T. F. & Houle, D. (2009) Corrigendum. J. Evol. Biol. 22:913-915.
#' @author Geir H. Bolstad
#' @examples
#' # Simulating a posterior distribution
#' # (or bootstrap distribution) of a G-matrix:
#' G <- matrix(c(1, 1, 0, 1, 4, 1, 0, 1, 2), ncol = 3)
#' G_mcmc <- sapply(c(G), function(x) rnorm(10, x, 0.01))
#' G_mcmc <- t(apply(G_mcmc, 1, function(x) {
#'   G <- matrix(x, ncol = sqrt(length(x)))
#'   G[lower.tri(G)] <- t(G)[lower.tri(G)]
#'   c(G)
#' }))
#'
#' # Simulating a posterior distribution
#' # (or bootstrap distribution) of trait means:
#' means <- c(1, 1.4, 2.1)
#' means_mcmc <- sapply(means, function(x) rnorm(10, x, 0.01))
#'
#' # Mean standardizing the G-matrix:
#' G_mcmc <- meanStdGMCMC(G_mcmc, means_mcmc)
#'
#' # Estimating average evolvability paramters:
#' evolvabilityMeansMCMC(G_mcmc)
#' @keywords array algebra
#' @importFrom stats median var
#' @export
evolvabilityMeansMCMC <- function(G_mcmc) {
  check_g_mcmc(G_mcmc)
  X <- list()
  results <- apply(G_mcmc, 1, function(G) {
    est <- mean_evolvability(matrix(G, ncol = sqrt(length(G))))
    c(est$estimates[c("e_mean", "e_min", "e_max", "r_mean", "c_mean", "a_mean",
                      "i_mean")],
      singular = est$singular, not_pd = est$not_pd)
  })
  results <- t(results)
  failed <- results[, "singular"] == 1
  warn_not_invertible(failed)
  warn_not_positive_definite(results[, "not_pd"] == 1)
  X$post.dist <- coda::as.mcmc(results[, seq_len(7), drop = FALSE])
  conditional <- c("c_mean", "a_mean", "i_mean")
  X$post.medians <- summarize_posterior_with_failed_draws(
    X$post.dist, invalid = any(failed) & colnames(X$post.dist) %in% conditional
  )
  X$call <- match.call()
  class(X) <- "evolvabilityMeansMCMC"
  X
}

#' @export
print.evolvabilityMeansMCMC <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nPosterior medians and 95% HPD intervals:\n")
  print(x$post.medians)
}
