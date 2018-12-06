#' Simulate responses from \code{\link{Almer}} fit
#' 
#' \code{simulate_Almer} Simulate responses from \code{\link{Almer}} fit
#' 
#' \code{simulate_Almer} ...
#' 
#' @param mod fitted object from \code{\link{Almer}}
#' @param nsim numder of simulations.
#'  
#' @details This function is only included as the simulate.merMod funciton did not seem to work properly when the 
#' number of random effect levels equal the number of observations. 
#' 
#' @return \code{simulate_Almer} a matrix of simulated responses, colunms correspond to each simulations.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @importFrom lme4 VarCorr ranef fixef 
#' @importFrom Matrix t
#' 
#' @export
#' 

# change name to Almer_sim
# check if this works with Almer_SE also (I think it should)
simulate_Almer <- function(mod, nsim = 10){
  samp_size_ranef <- c(sapply(ranef(mod), function(x) length(x[,1])))
  SD <- as.data.frame(VarCorr(mod))
  ranef_samples_list <- list()
  for(i in names(samp_size_ranef)){
    rn <- rnorm(samp_size_ranef[i]*nsim, 0 , SD[SD$grp == i,"sdcor"])
    ranef_samples_list[[i]] <- matrix(rn, ncol = nsim)
  }
  ranef_samples <- do.call(rbind, ranef_samples_list)
  Ranef <- do.call(cbind, apply(ranef_samples, 2, function(x) t(mod@pp$Zt)%*%x))
  Residuals <- rnorm(length(residuals(mod))*nsim, 0 , SD[SD$grp == "Residual","sdcor"])
  Fixef <- matrix(rep(mod@pp$X%*%fixef(mod), nsim), ncol = nsim)
  y <- Fixef + Ranef + Residuals
  colnames(y) <- paste("Sim", 1:nsim, sep = "_")
  return(y)
}


