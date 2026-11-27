#' Computing a conditional sub-matrix of G
#'
#' \code{conditinoalG} calculates a conditional variance matrix.
#'
#' @param G A symmetric variance matrix (see Details).
#' @param condition_on Either an integer with the column number indicating which
#'   trait to condition on or a vector with several column numbers (integers).
#' @details The function calculates a sub-matrix of \code{G} conditional on the
#'   traits defined by the the \code{condition_on} vector. The function is based
#'   on equation 3 in Hansen et al. (2003). An error is given if \code{G} is not
#'   symmetric. If the variance matrix of the traits in \code{condition_on}
#'   cannot be inverted by \code{solve()}, a warning is issued and a matrix of
#'   \code{NA} is returned. A warning is also issued if \code{G} is not positive
#'   definite (it has negative eigenvalues); the result is then still
#'   calculated but may not be meaningful.
#' @return A matrix that is a sub-matrix of the input matrix conditional on the
#'   non-included traits.
#' @references Hansen TF, Armbruster WS, Carlsson ML & Pélabon C. 2003.
#'   Evolvability and genetic constraint in Dalechampia blossoms: genetic
#'   correlations and conditional evolvability. J. Exp. Zool. 296B:23-39.
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
#' # Computing a conditional 2x2 sub-matrix by conditioning on
#' # trait 3 and 4:
#' G_sub_conditional <- conditionalG(G, condition_on = c(3, 4))
#' G_sub_conditional
#'
#' # The average evolvabilities of this matrix can then be
#' # compared can than be compared to the average evolvabilities
#' # of the corresponding unconditional sub-matrix of G:
#' evolvabilityMeans(G_sub_conditional)
#' evolvabilityMeans(G[-c(3, 4), -c(3, 4)])
#' @keywords array algebra
#' @export
conditionalG <- function(G, condition_on = NULL) {
  if (is.null(condition_on)) {
    stop("Specify one or several column numbers in a vector for 
          the argument condition_on (e.g. condition_on = c(1, 3) 
          for column 1 and 3)")
  }
  check_g(G)
  if (!is.numeric(condition_on) ||
      any(condition_on < 1 | condition_on > nrow(G) |
          condition_on != round(condition_on))) {
    stop("condition_on must contain column numbers of G.")
  }
  condition_on <- sort(unique(condition_on))
  if (length(condition_on) == nrow(G)) {
    stop("condition_on cannot include all the columns of G.")
  }
  var_y <- G[-condition_on, -condition_on, drop = FALSE]
  var_x <- G[condition_on, condition_on, drop = FALSE]
  cov_yx <- G[-condition_on, condition_on, drop = FALSE]
  warn_not_positive_definite(is_not_positive_definite(G))
  var_x_inv <- solve_g(var_x)
  if (is.null(var_x_inv)) {
    warn_not_invertible(
      TRUE, "The variance matrix of the traits in condition_on",
      "the conditional variance matrix is NA"
    )
    return(matrix(NA_real_, nrow(var_y), ncol(var_y),
                  dimnames = dimnames(var_y)))
  }
  return(var_y - cov_yx %*% var_x_inv %*% t(cov_yx))
}
