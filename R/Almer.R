#' Linear mixed model with correlated random effects structure
#' 
#' \code{Almer} Linear mixed model with correlated random effects structure
#' 
#' \code{Almer} fits a univariate linear mixed model incorporating a correlated random effects structure. Can be used to fit 
#' phylogenetic mixed models and animal models. The function is based on the \code{\link{lme4}} package and is very similar to \code{\link{lmer}}.
#' 
#' @param A an optional named list of sparce matrices. The names must correspond to the names of the random effects in the formula argument. All levels of the random effect should appear as row and column names for the matrices (SAME ORDER? or is names sufficient?)
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
#' @importFrom Matrix Matrix t
#' 
#' @export



Almer <- function(formula, data = NULL, REML = TRUE, A = list(), 
                  control = lmerControl(check.nobs.vs.nlev  = "ignore", 
                                        check.nobs.vs.rankZ = "ignore", 
                                        check.nobs.vs.nRE   = "ignore"),
                  start = NULL, verbose = 0L, subset, weights, na.action,
                  offset, contrasts = NULL, devFunOnly = FALSE, ...){
  cholA <- lapply(A, chol)
  mod <- lFormula(formula, data, REML, subset, weights, na.action, offset, contrasts, control, ...)
  for(i in seq_along(cholA)){
    j <- match(names(cholA)[i], names(mod$reTrms$cnms))
    if(length(j)>1) stop("an A matrix can only be associated with one random effect term")
    mod$reTrms$Ztlist[[j]] <- cholA[[i]] %*% mod$reTrms$Ztlist[[j]]
  }
  mod$reTrms$Zt <- do.call(rbind, mod$reTrms$Ztlist)
  devfun <- do.call(mkLmerDevfun, mod)
  opt <- optimizeLmer(devfun, optimizer = "Nelder_Mead", ...)
  mkMerMod(environment(devfun), opt = opt, reTrms = mod$reTrms, fr = mod$fr)
}