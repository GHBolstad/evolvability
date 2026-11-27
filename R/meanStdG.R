#' Mean standardize a variance matrix
#'
#' \code{meanStdG} mean standardizes a variance matrix (e.g. a G-matrix).
#'
#' @param G A symmetric variance matrix. An error is given if it is not
#'   symmetric.
#' @param means A vector of trait means.
#' @return A mean standardized variance matrix.
#' @details An error is given if \code{G} is not symmetric, or if the length of
#'   \code{means} differs from the dimension of \code{G}.
#' @author Geir H. Bolstad
#' @examples
#' G <- matrix(c(1, 1, 0, 1, 4, 1, 0, 1, 2), ncol = 3)
#' means <- c(1, 1.4, 2.1)
#' meanStdG(G, means)
#' @keywords array algebra
#' @export
meanStdG <- function(G, means) {
  standardize_g(G, means)
}
