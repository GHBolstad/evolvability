#' Calculate posterior distribution of evolvability parameters from a set of
#' selection gradients
#'
#' \code{evolvabilityBetaMCMC} calculates (unconditional) evolvability (e),
#' respondability (r), conditional evolvability (c), autonomy (a) and
#' integration (i) from selection gradients given the posterior distribution of
#' an additive-genetic variance matrix. These measures and their meanings are
#' described in Hansen and Houle (2008).
#'
#' @param G_mcmc posterior distribution of a variance matrix in the form of a
#'   table. Each row in the table must be one iteration of the posterior
#'   distribution (or bootstrap distribution). Each iteration of the matrix must
#'   be on the form as given by \code{c(x)}, where \code{x} is a matrix. A
#'   posterior distribution of a matrix in the slot \code{VCV} of a object of
#'   class \code{MCMCglmm} is by default on this form.
#' @param Beta either a vector or a matrix of unit length selection gradients
#'   stacked column wise. A warning is issued if any gradient's Euclidean
#'   length differs from one by more than \code{sqrt(.Machine$double.eps)} or
#'   is not finite. Gradients are not automatically normalized.
#' @param post.dist logical: should the posterior distribution of the
#'   evolvability parameters be saved.
#' @details An error is given if any draw in \code{G_mcmc} is not symmetric. If
#'   \code{solve(G)} fails for a posterior draw because the matrix is
#'   singular or computationally singular, a warning is issued. Conditional
#'   evolvability, autonomy and integration are \code{NA} for that draw, and
#'   their posterior summaries are \code{NA}; the other results are still
#'   calculated. A warning is also issued if draws are not positive definite
#'   (they have negative eigenvalues); the results are then still calculated
#'   but may not be meaningful for those draws.
#' @return An object of \code{class} \code{'evolvabilityBetaMCMC'}, which is a
#'   list with the following components:
#' \tabular{llllll}{
#' \code{eB} \tab\tab\tab\tab The posterior median and highest posterior density
#'  interval of evolvability for each selection gradient. \cr
#' \code{rB} \tab\tab\tab\tab The posterior median and highest posterior density
#'  interval of respondability for each selection gradient. \cr
#' \code{cB} \tab\tab\tab\tab The posterior median and highest posterior density
#'  interval of conditional evolvability for each selection gradient. \cr
#' \code{aB} \tab\tab\tab\tab The posterior median and highest posterior density
#'  interval of autonomy for each selection gradient.\cr
#' \code{iB} \tab\tab\tab\tab The posterior median and highest posterior density
#'  interval of integration for each selection gradient.\cr
#' \code{Beta} \tab\tab\tab\tab The matrix of selection gradients. \cr
#' \code{summary} \tab\tab\tab\tab The means of evolvability parameters across
#'  all selection gradients. \cr
#' \code{post.dist} \tab\tab\tab\tab The full posterior distribution.
#' }
#' @references Hansen, T. F. & Houle, D. (2008) Measuring and comparing
#'   evolvability and constraint in multivariate characters. J. Evol. Biol.
#'   21:1201-1219.
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
#' # Generating selection gradients in five random directions:
#' Beta <- randomBeta(5, 3)
#'
#' # Calculating evolvability parameters:
#' x <- evolvabilityBetaMCMC(G_mcmc, Beta, post.dist = TRUE)
#' summary(x)
#' @keywords array algebra multivariate
#' @importFrom stats median
#' @export
evolvabilityBetaMCMC <- function(G_mcmc, Beta, post.dist = FALSE) {
  check_beta_length(Beta)
  Beta <- cbind(Beta)
  if (nrow(Beta) != check_g_mcmc(G_mcmc)) {
    stop("Beta must have as many rows as the dimension of G.", call. = FALSE)
  }
  n <- ncol(Beta)
  X1 <- t(apply(G_mcmc, 1, function(G) {
    G <- matrix(G, ncol = sqrt(length(G)))
    p <- beta_parameters(G, Beta)
    c(eB = p$e, rB = p$r, cB = p$c, aB = p$a, iB = p$i,
      failed_solve = p$failed, not_pd = p$not_pd)
  }))
  failed <- X1[, 5 * n + 1] == 1
  not_pd <- X1[, 5 * n + 2] == 1
  X1 <- X1[, seq_len(5 * n), drop = FALSE]
  warn_not_invertible(failed)
  warn_not_positive_definite(not_pd)
  X2 <- list(eB = X1[, 1:n, drop = FALSE],
             rB = X1[, (n + 1):(2 * n), drop = FALSE],
             cB = X1[, (2 * n + 1):(3 * n), drop = FALSE],
             aB = X1[, (3 * n + 1):(4 * n), drop = FALSE],
             iB = X1[, (4 * n + 1):(5 * n), drop = FALSE]
             )
  X_summary <- cbind(sapply(X2, function(x) apply(x, 1, mean)))
  colnames(X_summary) <- c("e_mean", "r_mean", "c_mean", "a_mean", "i_mean")
  X <- lapply(seq_along(X2), function(i) {
    x <- X2[[i]]
    if (any(failed) && i >= 3) {
      t(summarize_posterior_with_failed_draws(x, rep(TRUE, ncol(x))))
    } else {
      rbind(median = apply(x, 2, median), t(coda::HPDinterval(coda::mcmc(x))))
    }
  })
  names(X) <- names(X2)
  X$Beta <- Beta
  X$summary <- if (any(failed)) {
    t(summarize_posterior_with_failed_draws(
      X_summary, c(FALSE, FALSE, TRUE, TRUE, TRUE)
    ))
  } else {
    rbind(median = apply(X_summary, 2, median),
          t(coda::HPDinterval(coda::mcmc(X_summary))))
  }
  X$call <- match.call()
  class(X) <- "evolvabilityBetaMCMC"
  X2$summary <- X_summary
  if (post.dist) X$post.dist <- X2
  return(X)
}

#' @export
print.evolvabilityBetaMCMC <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nEvolvability, posterior medians and 95% HPD intervals:\n")
  print(t(x$eB))
  cat("\nRespondability, posterior medians and 95% HPD intervals:\n")
  print(t(x$rB))
  cat("\nConditional evolvability, posterior medians and 95% HPD intervals:\n")
  print(t(x$cB))
  cat("\nAutonomy, posterior medians and 95% HPD intervals:\n")
  print(t(x$aB))
  cat("\nIntegration, posterior medians and 95% HPD intervals:\n")
  print(t(x$iB))
}


#' Summarizing posterior distribution of evolvability parameters over a set of
#' selection gradients
#'
#' \code{summary} method for class \code{'evolvabilityBetaMCMC'}.
#'
#' @param object an object of class \code{'evolvabilityBetaMCMC'}.
#' @param ... additional arguments affecting the summary produced.
#' @return A list with the following components:
#' \tabular{llllll}{
#'   \code{Averages} \tab\tab\tab\tab The averages of the evolvability parameters
#'    over all selection gradients. \cr
#'   \code{Minimum} \tab\tab\tab\tab The minimum (given by the posterior median)
#'    of the evolvability parameters over all selection gradients. \cr
#'   \code{Maximum} \tab\tab\tab\tab The maximum (given by the posterior median)
#'    of the evolvability parameters over all selection gradients.
#' }
#' @author Geir H. Bolstad
#' @seealso \code{\link{evolvabilityBetaMCMC}}
#' @keywords array algebra
#' @export
summary.evolvabilityBetaMCMC <- function(object, ...) {
  X <- list()
  X$call <- object$call
  X$Averages <- object$summary
  X$Minimum <- sapply(object[1:5], function(x) {
    if (anyNA(x[1, ])) rep(NA_real_, nrow(x)) else
      x[, which.min(x[1, ])]
  })
  colnames(X$Minimum) <- paste(colnames(X$Minimum), "_min", sep = "")
  X$Maximum <- sapply(object[1:5], function(x) {
    if (anyNA(x[1, ])) rep(NA_real_, nrow(x)) else
      x[, which.max(x[1, ])]
  })
  colnames(X$Maximum) <- paste(colnames(X$Maximum), "_max", sep = "")
  class(X) <- "summary.evolvabilityBetaMCMC"
  X
}

#' @export
print.summary.evolvabilityBetaMCMC <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nAverage values over all selection gradients, \n presented as posterior medians and 95% HPD intervals:\n")
  print(x$Averages)
  cat("\nMinimum (the direction with the lowest posterior median):\n")
  print(x$Minimum)
  cat("\nMaximum (the direction with the highest posterior median):\n")
  print(x$Maximum)
}
