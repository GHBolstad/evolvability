#' Linear mixed model for responsvariables with uncertainty
#' 
#' \code{Almer_SE} Linear mixed model for responsvariables with uncertainty
#' 
#' @param formula as in \code{\link{lmer}}.
#' @param SE a vector of standard errors associated with the response variable.
#' @param maxiter maximum number of iterations.
#' @param ... optional arguments, see \code{\link{Almer}}.
#' 
#' @return \code{Almer_SE} an object of class \code{\link{merMod}}.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @importFrom lme4 VarCorr lmerControl
#' 
#' @export

Almer_SE <- function(formula, SE=NULL, maxiter = 100, ...){
  
  if(is.null(SE)) stop("Use Almer instead")
  wgth <- 1/(1+SE^2)
  mod <- Almer(formula, weights = wgth, ...)

  for(i in 1:maxiter){
    rvariance <- attr(lme4::VarCorr(mod), "sc")^2
    wgth <- 1/(1+(SE^2)/rvariance)
    mod <- Almer(formula, weights = wgth, ...)
    if(abs(rvariance-attr(lme4::VarCorr(mod), "sc")^2)<1e-8) break()
  }
  if(i == maxiter){
    warning("Optimization of weights reached maximum number of iterations.")
  }
  
  return(mod)
}

