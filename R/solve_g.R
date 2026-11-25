solve_g <- function(G) {
  tryCatch(solve(G), error = function(e) {
    if (!grepl("system is (exactly|computationally) singular",
               conditionMessage(e))) {
      stop(e)
    }
    NULL
  })
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
