# Internal helpers shared by several functions. None of these are exported.

# Stops if G is not a (numerically) symmetric square matrix.
check_symmetric <- function(G, name = "G") {
  if (!isSymmetric(unname(as.matrix(G)), tol = sqrt(.Machine$double.eps))) {
    stop(name, " is not symmetric.", call. = FALSE)
  }
  invisible(TRUE)
}

# Checks that G is symmetric and mean standardizes it. The default
# `means = 1` means that G is left unstandardized.
standardize_g <- function(G, means = 1, name = "G") {
  check_symmetric(G, name)
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
# G cannot be inverted (see solve_g()).
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
       failed = is.null(G_inv))
}

# The average evolvability parameters of a G matrix (Hansen and Houle 2008),
# calculated from the eigenvalues of G.
mean_evolvability <- function(G) {
  values <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
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
  c(e_mean = e_mean, e_min = min(values), e_max = max(values),
    c_mean = c_mean, r_mean = r_mean, a_mean = a_mean, i_mean = 1 - a_mean)
}

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
