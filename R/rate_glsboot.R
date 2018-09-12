#' Boostrap of the rate gls model fit
#' 
#' \code{rate_glsboot} Boostrap of the rate gls model fit
#' 
#' \code{rate_glsboot} Provides a parametric boostrap of the generalized least squares rate model fit, using \code{\link{rate_sim}}  
#' 
#' @param mod output from \code{\link{rate_gls}}.
#' @param n number of bootsrap samples

rate_glsboot <- function(mod, n = 10){
  boot_distribution <- matrix(NA, ncol = 6, nrow = n)
  colnames(boot_distribution) <- c("Intercept", "Slope", "a", "b", "mean_x", "Sigma^2")
  
  for(i in 1:n){
    Data <- ratesim(tree=mod$tree, 
                    startv_x = mod$param["mean_x", ], 
                    sigma_x = sqrt(mod$param["Sigma^2", ]),
                    a = mod$param["a",],
                    b = mod$param["b",],
                    model = mod$model)
    boot_distribution[i,] <- c(rate_gls(x=Data$x, y=Data$y, tree=mod$tree, Beta = as.matrix(mod$param[1:2,]), model = mod$model, silent = TRUE)$param)
    print(paste("Bootstrap iteration", i))
  }
  
  return(cbind(mod$param, t(apply(boot_distribution, 2, function(x) quantile(x, probs = c(0.025, 0.975))))))
  
}

