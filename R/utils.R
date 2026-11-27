# Internal helpers shared by several functions. None of these are exported.

#### Checks of the input ####

# Warns if selection gradients are not of unit length (columns of Beta, or rows
# if `rows = TRUE`).
check_beta_length <- function(Beta, rows = FALSE) {
  Beta <- cbind(Beta)
  lengths <- if (rows) sqrt(rowSums(Beta^2)) else sqrt(colSums(Beta^2))
  if (any(!is.finite(lengths) | abs(lengths - 1) > sqrt(.Machine$double.eps))) {
    warning("Selection gradients (in Beta) are not of unit length", call. = FALSE)
  }
}

# Stops if G is not a numeric, finite and (numerically) symmetric square matrix.
check_g <- function(G, name = "G") {
  G <- as.matrix(G)
  if (!is.numeric(G) || any(!is.finite(G))) {
    stop(name, " must be numeric without missing or non-finite values.",
         call. = FALSE)
  }
  if (!isSymmetric(unname(G), tol = sqrt(.Machine$double.eps))) {
    stop(name, " is not symmetric.", call. = FALSE)
  }
  invisible(TRUE)
}

# Stops if the rows of G_mcmc (one vectorized G matrix per posterior draw) do
# not all describe numeric, finite and symmetric square matrices. Returns the
# dimension of the G matrices.
check_g_mcmc <- function(G_mcmc, name = "G_mcmc") {
  G_mcmc <- as.matrix(G_mcmc)
  k <- sqrt(ncol(G_mcmc))
  if (k != round(k)) {
    stop("The number of columns of ", name, " must be the square of the ",
         "dimension of G (one column per element of G).", call. = FALSE)
  }
  if (!is.numeric(G_mcmc) || any(!is.finite(G_mcmc))) {
    stop(name, " must be numeric without missing or non-finite values.",
         call. = FALSE)
  }
  transposed <- as.vector(t(matrix(seq_len(k^2), k)))
  asymmetric <- rowSums(abs(G_mcmc - G_mcmc[, transposed, drop = FALSE]) >
                          sqrt(.Machine$double.eps) * max(abs(G_mcmc))) > 0
  if (any(asymmetric)) {
    stop(name, " contains G matrices that are not symmetric (first in row ",
         which(asymmetric)[1], ").", call. = FALSE)
  }
  k
}

#### Inversion and positive definiteness of G ####

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

# TRUE if the eigenvalues show that G is not positive (semi-)definite, i.e.
# there are negative eigenvalues beyond numerical noise. Singular but positive
# semi-definite matrices are not flagged (see solve_g() for those).
has_negative_eigenvalues <- function(values) {
  min(values) < -sqrt(.Machine$double.eps) * max(abs(values))
}

is_not_positive_definite <- function(G) {
  has_negative_eigenvalues(
    eigen(G, symmetric = TRUE, only.values = TRUE)$values
  )
}

# Warns that G (or `sum(failed)` of the posterior draws of G) could not be
# inverted. `failed` is a logical vector with one element per G matrix: a
# single element for one matrix, one per draw in the MCMC functions.
warn_not_invertible <- function(failed, name = "G",
                                consequence = paste("conditional evolvability,",
                                                    "autonomy and integration",
                                                    "are NA")) {
  if (!any(failed)) return(invisible())
  if (length(failed) == 1) {
    warning(name, " could not be inverted by solve(); ", consequence, ".",
            call. = FALSE)
  } else {
    warning(sum(failed), " G matrix/matrices could not be inverted by solve(); ",
            consequence, " for those draws and their posterior summaries.",
            call. = FALSE)
  }
}

# Warns that G (or `sum(not_pd)` of the posterior draws of G) is not positive
# definite. Results are still calculated.
warn_not_positive_definite <- function(not_pd, name = "G") {
  if (!any(not_pd)) return(invisible())
  if (length(not_pd) == 1) {
    warning(name, " is not positive definite (it has negative eigenvalues); ",
            "the results may not be meaningful.", call. = FALSE)
  } else {
    warning(sum(not_pd), " G matrix/matrices are not positive definite (they ",
            "have negative eigenvalues); the results may not be meaningful ",
            "for those draws.", call. = FALSE)
  }
}

# Posterior medians and HPD intervals of the columns of x. Columns flagged in
# `invalid` (e.g. parameters that are NA for failed draws) are returned as NA.
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

#### Calculation of the evolvability parameters ####

# Checks that G is symmetric and mean standardizes it. The default
# `means = 1` means that G is left unstandardized.
standardize_g <- function(G, means = 1, name = "G") {
  check_g(G, name)
  k <- nrow(G)
  if (length(means) == 1 && means == 1) {
    means <- rep(1, k)
  }
  if (length(means) != k) {
    stop("The length of the means must equal the dimension of ", name, ".",
         call. = FALSE)
  }
  G / (means %*% t(means))
}

# Evolvability, respondability, conditional evolvability, autonomy and
# integration of the (column) selection gradients in Beta (a matrix), given a
# symmetric G. Conditional evolvability, autonomy and integration are NA when
# G cannot be inverted (see solve_g()); the other parameters are still
# calculated.
beta_parameters <- function(G, Beta) {
  GB <- G %*% Beta
  e_B <- colSums(Beta * GB)
  r_B <- sqrt(colSums(GB^2))
  G_inv <- solve_g(G)
  c_B <- if (is.null(G_inv)) {
    rep(NA_real_, ncol(Beta))
  } else {
    1 / colSums(Beta * (G_inv %*% Beta))
  }
  a_B <- c_B / e_B
  list(e = e_B, r = r_B, c = c_B, a = a_B, i = 1 - a_B,
       failed = is.null(G_inv), not_pd = is_not_positive_definite(G))
}

# The average evolvability parameters of a G matrix (Hansen and Houle 2008),
# calculated from the eigenvalues of G. Returns a list with the estimates, and
# whether G is singular (see solve_g()) or not positive definite. The
# conditional evolvability, autonomy and integration are NA if G is singular;
# the other parameters are still calculated.
mean_evolvability <- function(G) {
  values <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  singular <- rcond(as.matrix(G)) < .Machine$double.eps
  VAR <- function(x) sum((x - mean(x))^2) / length(x)
  k <- nrow(G)
  e_mean <- mean(values)
  Heig <- 1 / mean(1 / values)
  Ieig <- VAR(values) / e_mean^2
  Ieig2 <- VAR(values^2) / mean(values^2)^2
  Iinveig <- VAR(1 / values) / mean(1 / values)^2
  c_mean <- Heig * (1 + (2 * Iinveig) / (k + 2))
  r_mean <- sqrt(mean(values^2)) * (1 - (Ieig2 / (4 * (k + 2))))
  a_mean <- (Heig / e_mean) *
    (1 + 2 * (Ieig + Iinveig - 1 + (Heig / e_mean) +
                2 * Ieig * Iinveig / (k + 2)) / (k + 2))
  out <- c(e_mean = e_mean, e_min = min(values), e_max = max(values),
           c_mean = c_mean, r_mean = r_mean, a_mean = a_mean,
           i_mean = 1 - a_mean)
  if (singular) {
    out[c("c_mean", "a_mean", "i_mean")] <- NA_real_
  }
  list(estimates = out, singular = singular,
       not_pd = has_negative_eigenvalues(values))
}

#### Helpers for the rate models and the mixed models ####

# Replaces eigenvalues below `tol` with `tol` to ensure a positive definite
# matrix.
make_positive_definite <- function(M, tol = 1e-8) {
  e <- eigen(M)
  if (any(e$values < tol)) {
    e$values[e$values < tol] <- tol
    M <- e$vectors %*% diag(e$values) %*% t(e$vectors)
  }
  M
}

# Name of the temporary response column used when refitting a model to
# simulated data (see refit_simulated()).
sim_response_name <- "x86eb0d23_ebe6_4b37_bbc7_1d164b2f7c26"

# Simulates `nsim` new responses from the fitted Almer model `mod` and refits
# the model to each of them. Returns a list of fitted models.
refit_simulated <- function(mod, nsim) {
  y <- Almer_sim(mod, nsim)
  dt <- mod@frame
  new_formula <- stats::as.formula(paste(sim_response_name, "~ ."))
  fits <- lapply(seq_len(ncol(y)), function(i) {
    dt[[sim_response_name]] <- y[, i]
    update(mod, new_formula, data = dt)
  })
  names(fits) <- colnames(y)
  fits
}
