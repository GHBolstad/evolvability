#' Computing a conditional sub matrix of G
#'
#' \code{conditinoalG} calculates a conditional variance matrix.
#'
#' @param G A variance matrix (must be symmetric and positive definite).
#' @param condition_on Either an integer with the column number on which to
#'   condition on or a vector with several column numbers (integers).
#' @details The function is motivated by the conditional evolvability, but in a
#'   multivariate setting. It can be used on any symmetric positive definite
#'   variance matrix. In the computation, the matrix is first rotated so that
#'   the all traits, but the one specified by the `condition_on` vector, are
#'   uncorrelated. Second, the \code{evolvabilityBeta} function is used to
#'   extract the conditional evolvabilities of the diagonal. Third, the traits
#'   specified by the `condition_on` vector are removed and the resulting
#'   sub-matrix is rotated back to reflect the original traits in `G`. The
#'   returned variance matrix is a sub-matrix of `G` that contain the variance
#'   that is independent of the traits specified by `condition_on`.
#' @return A matrix that is a sub-matrix of G conditional on the non-included traits.
#' @author Geir H. Bolstad
#' @examples
#' # Constructing a G-matrix:
#' G <- matrix(c(
#'   1, 1, 0, 1,
#'   1, 2, 1, 1,
#'   0, 1, 2, 1,
#'   1, 1, 1, 3
#' ), ncol = 4)
#'
#' # Computing a conditional 2x2 sub-matrix by conditioning on trait 1 and 3:
#' G_sub_conditional <- conditionalG(G, condition_on = c(1, 3))
#'
#' # The average evolvabilities of this matrix can then be compared can than be
#' # compared to the average evolvabilities of the corresponding unconditional
#' # sub matrix of G:
#' evolvabilityMeans(G_sub_conditional)
#' evolvabilityMeans(G[-c(1, 3), -c(1, 3)])
#' @keywords array algebra
#' @export
conditionalG <- function(G, condition_on = NULL) {
  if (is.null(condition_on)) {
    stop("Specify on or several column numbers in a vector (e.g. c(1,3) for column 
         1 and 3)")
  }
  if (any(G[upper.tri(G)] != G[t(lower.tri(G))])) {
    stop("G is not symmetric.")
  }
  n <- ncol(G)
  z_sub <- c(1:n)[-condition_on]
  G_sub <- G[-condition_on, -condition_on]
  eVec_sub <- eigen(G_sub)$vectors
  eVec <- diag(n)
  eVec[z_sub, z_sub] <- eVec_sub
  G_rotated <- t(eVec) %*% G %*% eVec
  conditionalG_rotated <- diag(evolvabilityBeta(G_rotated, diag(n))$c[z_sub])
  conditionalG <- eVec_sub %*% conditionalG_rotated %*% t(eVec_sub)
  return(conditionalG)
}
