#' Phylogenetic heritability
#'
#' \code{phylH} calculates the phylogenetic heritability from an \code{Almer}
#' model fit and provides associated uncertainty using parametric
#' bootstrapping.
#'
#' @param mod An object of class \code{'merMod'}
#' @param numerator The name of phylogenetic effect level
#' @param residual name of the residual effect level
#' @param nsim number of bootstraps
#' @return \code{phylH} returns a list with the REML estimate, the 95\%
#'   confidence interval from the parametric bootstrap, and the bootstrap
#'   samples.
#' @author Geir H. Bolstad
#' @examples
#' # See the vignette 'Phylogenetic mixed model'.
#' @importFrom lme4 VarCorr
#'
#' @export

phylH <- function(mod, numerator, residual = "Residual", nsim = 10) {
  boot <- vapply(refit_simulated(mod, nsim),
                 function(m) H(m, numerator, residual), numeric(1),
                 USE.NAMES = TRUE)
  return(list(phylH = c(H(mod, numerator, residual), 
                        quantile(boot, c(0.025, 0.975))),
              bootstrap = boot
              )
        )
}

H <- function(mod, numerator, residual = "Residual") {
  num <- attr(VarCorr(mod)[[numerator]], "stddev")^2
  if (residual == "Residual") {
    denom <- num + attr(VarCorr(mod), "sc")^2
  } else {
    denom <- num + attr(VarCorr(mod)[[residual]], "stddev")^2
  }
  H <- num / denom
  names(H) <- "Phylo Heritability"
  return(H)
}
