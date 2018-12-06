#' Phylogenetic heritability
#' 
#' \code{phylH} Phylogenetic heritability from a Almer fit
#' 
#' \code{phylH} calculates the phylogenetic heritability with uncertainty using parameteric bootstrapping
#' 
#' @param mod a merMod object
#' @param numerator name of phylogenetic effect level
#' @param residual name of the residual effect level  
#' @param nsim number of bootstraps
#' 
#' @return \code{phylH} returns a list with the REML estimate, the 95% confidence interval from the parametric bootstrap, and the bootstrap samples.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @importFrom lme4 VarCorr 
#' 
#' @export


phylH <- function(mod, numerator, residual = "Residual", nsim = 10){
  y <- simulate_Almer(mod, nsim)
  dt <- mod@frame
  boot <- apply(y, 2, function(x){
    dt$hopefullynotavariablealready2234uoiu234 <- x
    mod_updt <- update(mod, hopefullynotavariablealready2234uoiu234 ~ ., data = dt)
    H(mod_updt, numerator, residual)
  })
  return(list(phylH = c(H(mod, numerator, residual), quantile(boot, c(0.025, 0.975))),
         bootstrap = boot))
}


H <- function(mod, numerator, residual = "Residual"){
  num <- attr(VarCorr(mod)[[numerator]], "stddev")^2
  if(residual == "Residual") denom <- num + attr(VarCorr(mod), "sc")^2
  else denom <- num + attr(VarCorr(mod)[[residual]], "stddev")^2
  H <- num/denom
  names(H) <- "Phylo Heritability"
  return(H)
}

