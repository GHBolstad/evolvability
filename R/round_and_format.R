#' Rounds and formats in the same function
#'
#' \code{round_and_format} rounds and formats a numeric vector. This is useful for
#'  providing output for tables or plots in a standardized format.
#'
#' @param x A numeric vector.
#' @param digits Number of decimal places.
#' @param sign_digits Number of significant digits (if given this overrides
#'  \code{digits}).
#' @param scientific logical: whether encoding should be in scientific notation or
#'  not.
#' @param trim logical: if leading blanks for justification to common width should
#'  be excluded or not.
#' @return Rounded and formatted values as characters.
#' @author Geir H. Bolstad
#' @export
round_and_format <- function(x, digits = 2, sign_digits = NULL, scientific = FALSE,
                             trim = TRUE) {
  if (is.null(sign_digits)) {
    format(round(x, digits = digits),
      nsmall = digits, scientific = scientific,
      trim = trim
    )
  } else {
    x <- signif(x, digits = sign_digits)
    # Number of decimals needed to display sign_digits significant digits
    nsmall <- vapply(x, function(xi) {
      if (!is.finite(xi) || xi == 0) return(0)
      min(20, max(0, sign_digits - 1 - floor(log10(abs(xi)))))
    }, numeric(1))
    out <- vapply(seq_along(x), function(i) {
      format(x[i], nsmall = nsmall[i], scientific = scientific, trim = trim)
    }, character(1))
    names(out) <- names(x)
    out
  }
}
