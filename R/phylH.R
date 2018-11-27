#' Phylogenetic heritability
#' 
#' \code{phylH_boot} Phylogenetic heritability from a merMod object
#' 
#' \code{phylH_boot} calculates the phylogenetic heritability with uncertainty using parameteric bootstrapping
#' 
#' @param mod a mer mode object
#' @param numerator name of phylogenetic effect level
#' @param residual name of the residual effect level  
#' @param nsim number of bootstraps
#' 
#' @return \code{phylH_boot} 
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @importFrom lme4 VarCorr
#' 
#' @export



phylH_boot <- function(mod, numerator, residual = "Residual", nsim = 999){
  phylH(mod, numerator, residual)
  y <- simulate(mod, nsim)
  boot_samples <- NULL
  for(i in 1:nsim){
    mod_updt <- update(mod, y[,i] ~ .)
    boot_samples[i] <- phylH(mod_updt, numerator, residual)
  }
  return(list(phylH = c(phylH(mod, numerator, residual), quantile(boot_samples, c(0.025, 0.975)))),
              bootstrap_samples = boot_samples)
  
}


phylH <- function(mod, numerator, residual = "Residual"){
  num <- attr(VarCorr(mod)[[numerator]], "stddev")^2
  if(residual == "Residual") denom <- num + attr(VarCorr(mod), "sc")^2
  else denom <- num + attr(VarCorr(mod)[[residual]], "stddev")^2
  H <- num/denom
  names(H) <- "Phylo Heritability"
  return(H)
  }

