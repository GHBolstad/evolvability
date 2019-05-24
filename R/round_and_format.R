#' Rounds and formats in the same function
#'
#' \code{round_and_format} rounds and formats data
#'
#' @param x data in numeric format
#' @param digits Number of decimal places
#' @param sign_digits Number of significant digits (if given this overrides \code{digits})
#' @param scientfic Logical specifying whether encoding in scientific notation
#' @return rounded and formatted values
#' @author Geir H. Bolstad
#'
#' @export

round_and_format <- function(x, digits = 2, sign_digits = NULL, scientific = FALSE, trim = TRUE){
  if(is.null(sign_digits)) format(round(x, digits = digits), digits = digits, nsmall = digits, scientific = scientific, trim = trim)
  else format(signif(x, digits = sign_digits), digits = sign_digits, nsmall = sign_digits, scientific = scientific, trim = trim)
}