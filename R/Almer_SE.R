#' !!! NEED TO FIX !!! Linear mixed model with correlated random effects structure
#' 
#' \code{Almer_SE} !!! NEED TO FIX !!! Linear mixed model with correlated random effects structure
#' 
#' @param SE a vector of standard errors associated with the response variable.
#' @param A as in \code{\link{Almer}}.
#' @param formula as in \code{\link{lmer}}.
#' @param data as in \code{\link{lmer}}.
#' @param REML as in \code{\link{lmer}}.
#' @param control as in \code{\link{lmer}}.
#' @param start as in \code{\link{lmer}}.
#' @param verbose as in \code{\link{lmer}}.
#' @param subset as in \code{\link{lmer}}.
#' @param weights as in \code{\link{lmer}}.
#' @param na.action as in \code{\link{lmer}}.
#' @param offset as in \code{\link{lmer}}.
#' @param contrasts as in \code{\link{lmer}}.
#' @param devFunOnly as in \code{\link{lmer}}.
#' @param ... as in \code{\link{lmer}}.
#'  
#' 
#' @return \code{Almer_SE} an object of class \code{\link{merMod}}.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @importFrom lme4 VarCorr
#' 
#' @export

#### need to check if this is working ###

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
  return(mod)
}

