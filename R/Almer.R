#' Linear mixed model with correlated random effects structure
#' 
#' \code{Almer} Linear mixed model with correlated random effects structure
#' 
#' \code{Almer} fits a univariate linear mixed model incorporating a correlated random effects structure. Can be used to fit 
#' phylogenetic mixed models and animal models. The function is based on the \code{\link{lme4}} package and is very similar to \code{\link{lmer}},
#' apart from the A argument.
#' 
#' @param A an optional named list of sparce matrices. The names must correspond to the names of the random effects in the formula argument. 
#' All levels of the random effect should appear as row and column names for the matrices.
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
#' @return \code{Almer} an object of class \code{\link{merMod}}.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @importFrom lme4 lFormula lmerControl mkLmerDevfun optimizeLmer mkMerMod
#' @importFrom Matrix Matrix t chol
#' 
#' @export

Almer <- function(formula, data = NULL, REML = TRUE, A = list(), 
                  control = lmerControl(check.nobs.vs.nlev  = "ignore", 
                                        check.nobs.vs.rankZ = "ignore", 
                                        check.nobs.vs.nRE   = "ignore"),
                  start = NULL, verbose = 0L, subset, weights, na.action,
                  offset, contrasts = NULL, devFunOnly = FALSE, ...){
  mc <- match.call()
  cholA <- lapply(A, chol)
  mod <- lFormula(formula, data, control=control, ...) # have to check how to include the additional arguments
  for(i in seq_along(cholA)){
    j <- match(names(cholA)[i], names(mod$reTrms$cnms))
    if(length(j)>1) stop("an A matrix can only be associated with one random effect term")
    ranef_order <- match(rownames(mod$reTrms$Ztlist[[j]]), rownames(cholA[[i]]))
    mod$reTrms$Ztlist[[j]] <- cholA[[i]][ranef_order,ranef_order] %*% mod$reTrms$Ztlist[[j]]
  }
  mod$reTrms$Zt <- do.call(rbind, mod$reTrms$Ztlist)
  devfun <- do.call(mkLmerDevfun, mod)
  opt <- optimizeLmer(devfun, optimizer = "Nelder_Mead", ...)
  mod_obj <- mkMerMod(environment(devfun), opt = opt, reTrms = mod$reTrms, fr = mod$fr, mc)
  mod_obj@call <- evalq(mc)
  return(mod_obj)
}

