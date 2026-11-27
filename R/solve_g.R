# Inverts G, or returns NULL if G is (numerically) singular. The criterion is
# the same as in base::solve() (reciprocal condition number below machine
# epsilon), but it does not depend on the language of error messages. Other
# problems (e.g. non-square or non-finite G) still raise an error.
solve_g <- function(G) {
  G <- as.matrix(G)
  if (rcond(G) < .Machine$double.eps) {
    return(NULL)
  }
  solve(G)
}

summarize_posterior_with_failed_draws <- function(x, invalid) {
  result <- matrix(NA_real_, nrow = ncol(x), ncol = 3,
                   dimnames = list(colnames(x), c("median", "lower", "upper")))
  if (any(!invalid)) {
    valid <- x[, !invalid, drop = FALSE]
    result[!invalid, ] <- cbind(
      median = apply(valid, 2, median),
      coda::HPDinterval(coda::mcmc(valid))
    )
  }
  result
}
