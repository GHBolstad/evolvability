check_beta_length <- function(Beta, rows = FALSE) {
  Beta <- cbind(Beta)
  lengths <- if (rows) sqrt(rowSums(Beta^2)) else sqrt(colSums(Beta^2))
  if (any(!is.finite(lengths) | abs(lengths - 1) > sqrt(.Machine$double.eps))) {
    warning("Selection gradients (in Beta) are not of unit length", call. = FALSE)
  }
}
