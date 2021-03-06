% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/round_and_format.R
\name{round_and_format}
\alias{round_and_format}
\title{Rounds and formats in the same function}
\usage{
round_and_format(
  x,
  digits = 2,
  sign_digits = NULL,
  scientific = FALSE,
  trim = TRUE
)
}
\arguments{
\item{x}{A numeric vector.}

\item{digits}{Number of decimal places.}

\item{sign_digits}{Number of significant digits (if given this overrides
\code{digits}).}

\item{scientific}{logical: whether encoding should be in scientific notation or
not.}

\item{trim}{logical: if leading blanks for justification to common width should
be excluded or not.}
}
\value{
Rounded and formatted values as characters.
}
\description{
\code{round_and_format} rounds and formats a numeric vector. This is useful for
 providing output for tables or plots in a standardized format.
}
\author{
Geir H. Bolstad
}
