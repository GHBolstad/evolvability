#' Rounds and formats in the same function
#'
#' \code{round_and_format} rounds and formats data
#'
#' @param x data in numeric format
#' @param digits Number of decimals
#' @param scientfic Logical specifying whether encoding in scientific notation
#' @return rounded and formatted values
#' @author Geir H. Bolstad
#'
#' @export

round_and_format <- function(x, digits = 2, scientific = FALSE, trim = TRUE){
  format(round(x, digits = digits), digits = digits, nsmall = digits, scientific = scientific, trim = trim)
}