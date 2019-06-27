#' Generalized least squares rate model
#' 
#' \code{rate_gls} Generalized least squares rate model
#' 
#' \code{rate_gls} Fits a generalized least squares model to estimate parameters of 
#' the evolutoinary model of two traits x and y, where the evoluitonary rate of y 
#' depends on the value of x. Three models are implmented. In the two first, 
#' "predictor_BM" and "predictor_gBM", the evolution of y follows a Brownian motion 
#' with variance linear in x, while the evolution of x either follows a Brownian 
#' motion or a geometric brownian motion, respectively. In the third model, the 
#' residuals of the macroveoluitonary predictions of y have variance linear in x. 
#'
#' @param x explanatory variable must be equal to the length of y and tips on the 
#' tree. Note that the function mean centres x in the "predictor_BM" and "recent_evol" 
#' analyses, while it is mean stadardized (divided by the mean) in the "predictor_gBM".
#' @param y trait values of response variable.
#' @param species names of the species, must be equal in length and in the same order 
#' as x and y.
#' @param tree object of class \code{\link{phylo}}, needs to be ultrametric and with 
#' total length of unit, tips must have the same names as in \code{species}.
#' @param model either "predictor_BM", "predictor_gBM" or "recent_evol" (see Details).
#' @param startv vector of optional starting values for the a and b parameters.
#' (Ignored for the "recent_evol" model.)
#' @param maxiter maximum number of iterations for updating the GLS.
#' @param silent if the function should print the generalized sum of squares for each 
#' iteration.
#' @param useLFO acronym for use leave focal out when calulating the mean vector of 
#' the traits. In the "recent_evol" analysis, should the focal species be left out 
#' when calculating the corresponding species' mean. The correct way is to use TRUE, 
#' but in practice it has little effect and FALSE will speed up the model fit 
#' (particularely useful when bootstrapping) 
#' 
#' @details the generalized least squares (GLS) model fit a regression where 
#' explanatory variable is x (mean centred) and the response variable is vector of 
#' squared y values for the "predictor_BM" and "predictor_gBM" models and squared 
#' deviation from the evolutionary predictions (squared residuals) in the 
#' "recent_evol" model. From the intercept and slope of the GLS fit the following 
#' evolutionary paramters (a and b) are estimated. The evolutionary models are
#' described in a vignette and in Hansen et al. (in prep)
#' 
#' @return \code{rate_gls} returns a list with elements
#' \itemize{
#'  \item{"model"}{the fitted model "predictor_BM", "predictor_gBM" or 
#'  "recent_evolution}
#'  \item{"param"}{parameter estimates and standard errors, where "A" and 
#'  "B" are the intercept and slope of the GLS regression, "a" and "b" are parameters 
#'  of the evolutionary models (see vignette), and "sigma2_x" which is the BM-rate 
#'  for the "predictor_BM" and "predictor_gBM" models and the variance of x for the 
#'  "recent_evolution" model.}
#'  }
#'  \item{Rsquared}{The generalized r-square of the GLS model fit}
#'  \item{GLS_objective_score}{The score of the GLS objective function}
#'  \item{R}{Residual variance matrix}
#'  \item{tree}{The phylogenetic tree}
#'  \item{data}{The data used for the GLS regresssion}
#'  \item{convergence}{Whether the algoritm converged or not}
#'  \item{call}
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @importFrom ape vcv
#' @importFrom lme4 VarCorr
#' @importFrom Matrix Matrix
#' 
#' @export

# TODO: print.rate_gls function
# 

rate_gls <- function(x, y, species, tree, model = "predictor_BM", 
                     startv = list(a = NULL, b = NULL), maxiter = 100, silent = FALSE,
                     useLFO = TRUE){
  

  #### Phylogenetic relatedness matrix ####
  if(!ape::is.ultrametric(tree)) stop("The tree is not ultrametric") 
  A <- ape::vcv(tree)
  if(round(A[1,1], 5) != 1) stop("The tree is not standardized to unit length")
  A <- A[match(species, colnames(A)), match(species, colnames(A))] #ordering A according to species
  Ainv <- solve(A)

  #### y-variable ####
  # mean (weighted assuming BM process and the phylogeny)
  X <- matrix(rep(1, length(x)), ncol = 1) # design matrix
  mean_y <- solve(t(X)%*%Ainv%*%X)%*%t(X)%*%Ainv%*%y

  # mean centring
  y <- y - c(mean_y) # note that for the residual rate model, this does not have any effect as the squared deviation from the predictons are used

  #### x-variable ####
  # mean of x
  if(model == "predictor_BM")  mean_x <- solve(t(X)%*%Ainv%*%X)%*%t(X)%*%Ainv%*%x            
  if(model == "predictor_gBM") mean_x <- exp(solve(t(X)%*%Ainv%*%X)%*%t(X)%*%Ainv%*%log(x))      
  if(model == "recent_evol")   mean_x <- mean(x)

  # mean centring x or standardizing x
  if(model == "predictor_BM")  x <- x - c(mean_x)
  if(model == "predictor_gBM") x <- x/c(mean_x)
  if(model == "recent_evol")   x <- x - c(mean_x) 
    
  # Variance of x
  C <- t(chol(A)) #left cholesky factor
  if(model == "predictor_BM")  sigma2_x <- var(solve(C)%*%x) 
  if(model == "predictor_gBM") sigma2_x <- var(solve(C)%*%(log(x)-mean(log(x))))
  if(model == "recent_evol")   sigma2_x <- var(x)
  sigma2_x <- c(sigma2_x)
  sigma2_x_SE <- sqrt(2*sigma2_x^2/(length(x)+2)) # from Lynch and Walsh 1998 eq A1.10c

  #### Internal functions ####
  if(model == "predictor_BM"){
    a_func <- function(Beta, sigma2_x) Beta[1,1]
    b_func <- function(Beta, sigma2_x) 2*Beta[2,1]
    R_func <- function(a, b, sigma2_x, t_a, V=NULL) 2*a^2*t_a^2 + c(sigma2_x*(b/2)^2) * t_a*(1 - 2*t_a + 4*t_a^2)
    a_SE_func <- function(Beta_vcov, sigma2_x) sqrt(Beta_vcov[1,1])
    b_SE_func <- function(Beta_vcov, sigma2_x) sqrt(4*Beta_vcov[2,2])
  }
  if(model == "predictor_gBM"){
    f <- function(sigma2_x, t_a = 1){
      X <- try((exp(sigma2_x*t_a)+2*exp(-0.5*sigma2_x*t_a)-3)/(exp(sigma2_x*t_a)-1))
      X[is.na(X)] <- 0
      return(X)
    }
    a_func <- function(Beta, sigma2_x) Beta[1,1]-3*Beta[2,1]*f(sigma2_x)*(exp(sigma2_x/2)-1)
    b_func <- function(Beta, sigma2_x) 3*sigma2_x*Beta[2,1]/(2*f(sigma2_x))
    R_func <- function(a, b, sigma2_x, t_a, V=NULL){
      2*a^2*t_a^2 + 8*a*b*t_a/sigma2_x*(exp(0.5*sigma2_x*t_a)-1)+
        2*b^2/sigma2_x^2*(exp(sigma2_x*t_a)-1)*(exp(sigma2_x*t_a)-1+2*(exp(0.5*sigma2_x*t_a)-exp(0.5*sigma2_x))^2+
                                              4/3*(exp(0.5*sigma2_x)-exp(0.5*sigma2_x*t_a))*(f(sigma2_x, t_a)*exp(0.5*sigma2_x*t_a)-f(sigma2_x)*exp(0.5*sigma2_x))-
                                              4/9*f(sigma2_x, t_a)*exp(0.5*sigma2_x*t_a)*f(sigma2_x)*exp(0.5*sigma2_x)+
                                              2/9*f(sigma2_x)^2*exp(sigma2_x)
        )
    }
    a_SE_func <- function(Beta_vcov, sigma2_x) sqrt(Beta_vcov[1,1] + 9*f(sigma2_x)^2*(exp(sigma2_x/2)-1)^2*Beta_vcov[2,2] - 3*f(sigma2_x)*(exp(sigma2_x/2-1)*Beta_vcov[1,2]))
    b_SE_func <- function(Beta_vcov, sigma2_x) sqrt(9*sigma2_x^2*Beta_vcov[2,2]/(4*f(sigma2_x)^2))
  }
  if(model == "recent_evol"){
    a_func <- function(Beta, sigma2_x) Beta[1,1]  
    b_func <- function(Beta, sigma2_x) Beta[2,1]
    R_func <- function(a = NULL, b = NULL, sigma2_x, t_a, V){
      Vinv <- solve(V)
      dVmin1 <- solve(diag(x = diag(Vinv)))
      dVmin2 <- diag(x = diag(dVmin1)^2) #eqivalent to dVmin1%*%dVmin1 for diagonal matrces
      return(2* (dVmin2%*%Vinv) * (dVmin2%*%Vinv))
    }
    a_SE_func <- function(Beta_vcov, sigma2_x) sqrt(Beta_vcov[1,1]) 
    b_SE_func <- function(Beta_vcov, sigma2_x) sqrt(Beta_vcov[2,2]) 
  }

  
  
                  #######################
                  #### The GLS model ####
                  #######################
                  
  #### response variabe and initial values ####
  obj <- NULL                                     # objective function scores
  t_a <- as.vector(A[!lower.tri(A)])              # age of common ancestor        
  a   <- startv$a                                 # parameter of the process
  b   <- startv$b                                 # parameter of the process

  if(model != "recent_evol"){
    # response variable
    y2 <- y^2
    # finding starting values when not specified
    if(is.null(a) | is.null(b)){ 
      coef_lm <- coef(lm(y2~x))
      if(is.null(a)) a <- coef_lm[1]
      if(is.null(b)) b <- coef_lm[2]
    }
    R <- A                                          
    R[!lower.tri(R)] <- R_func(a, b, sigma2_x, t_a, V)
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
  }
  if(model == "recent_evol") {
    mod_Almer <- Almer(y ~ 1 + (1|species), A = list(species = Matrix::Matrix(A, sparse = TRUE)))
    sigma2_y <- VarCorr(mod_Almer)$species[1]
    residvar <- attr(VarCorr(mod_Almer), "sc")^2
    if(is.null(a)) a <- residvar 
    if(is.null(b)) b <- 0
    E <- c(a)*diag(length(x)) + diag(x = c(b*x))
    V <- sigma2_y*A + E
    R <- R_func(a, b, sigma2_x, t_a, V)
  }
  
  # inverse of residual covariance matrix
  Rinv <- solve(R)
  
  # gls design matrix  
  X <- as.matrix(cbind(rep(1, length(x)), x))
  
  #### iterative GLS ####
  for(i in 1:maxiter){
    if(model == "recent_evol") {
      y_predicted <- macro_pred(y=y, V=V, useLFO=useLFO)
      y2 <- (y-y_predicted)^2
    }
    
    # gls estimates
    Beta <- solve(t(X)%*%Rinv%*%X)%*%t(X)%*%Rinv%*%y2
    a <- a_func(Beta, sigma2_x)
    b <- b_func(Beta, sigma2_x)
    
    # updating residual covariance matrix using gls estimates
    if(model != "recent_evol") {
      R[!lower.tri(R)] <- R_func(a, b, sigma2_x, t_a)
      R[lower.tri(R)] <- t(R)[lower.tri(R)]
    }
    if(model == "recent_evol") {
      E <- c(a)*diag(length(x)) + diag(x = c(b*x))
      V <- sigma2_y*A + E
      R <- R_func(a, b, sigma2_x, t_a, V)
    }
    Rinv <- solve(R)
    
    # Objective function
    obj[i] <- t(y2-X%*%Beta)%*%Rinv%*%(y2-X%*%Beta)
    if(!silent) print(paste("Generalized sum of squares:", obj[i]))
    if(i>1){
      if(obj[i]<=obj[i-1] & obj[i-1]-obj[i]<1e-8) break()
    }
  }
  
  Beta_vcov <- solve(t(X)%*%Rinv%*%X)
  Beta_SE <- sqrt(diag(Beta_vcov))
  param <- cbind(rbind(Beta, a, b, sigma2_x), rbind(cbind(Beta_SE), a_SE_func(Beta_vcov, sigma2_x), b_SE_func(Beta_vcov, sigma2_x), sigma2_x_SE))
  colnames(param) <- c("Estimate", "SE")
  rownames(param) <- c("Intercept", "Slope", "a", "b", "sigma2_x")

  # R squared
  G <- matrix(rep(1, length(y2)), ncol = 1) # design matrix
  y_mean <- solve(t(G)%*%Rinv%*%G)%*%t(G)%*%Rinv%*%y2
  TSS <- t(y2-G%*%y_mean)%*%Rinv%*%(y2-G%*%y_mean)
  Rsquared <- 1-obj[length(obj)]/TSS

  if(i == maxiter){
    warning("Model reached maximum number of iterations without convergence")
    convergence <- "No convergence"
  }
  else  convergence <- "Convergence"
  
  report <-  list(model = model, param = param, Rsquared = Rsquared, GLS_objective_scores = obj,
                  R=R, tree = tree, data = list(y2 = y2, x = x, y = y), convergence = convergence)
  if(model == "recent_evol") report$sigma2_y <- sigma2_y
  class(report) = "rate_gls"
  report$call <- match.call()
  report
}


#' Plot of gls_rate object
#' 
#' \code{plot} Plots the gls rate regression 
#' 
#' \code{plot} ... 
#'
#' @param object a gls_rate object
#' @param scale either the variance scale ("VAR") or the standard deviation scale ("SD").
#' @param print_param logical if parameter estimates should be printed.
#' @param digits_param number of sigificant digits.
#' @param digits_rsquared number of decimal places
#' @param main as in \code{\link{plot}}.
#' @param xlab as in \code{\link{plot}}.
#' @param ylab as in \code{\link{plot}}.
#' @param col as in \code{\link{plot}}.
#' @param ... additional arguments passed to \code{\link{plot}}.
#' 
#' @details Plots the gls rate regression from a gls_rate object obtained from the gls_rate function. The regression line 
#' gives the expected variance or standard deviation (depending on scale).The regression is linear on the variance scale.
#'  
#' @return \code{plot} returns a plot of the gls rate regression
#' 
#' @author Geir H. Bolstad
#' 
#' @export

plot.rate_gls = function(object, scale = "SD", print_param = TRUE, digits_param = 2, digits_rsquared = 1, main = "GLS regression", xlab = "x", ylab = "Response", col = "grey",  ...){
  x <- seq(min(object$data$x), max(object$data$x), length.out = 100)
  y <- object$param["Intercept",1] + object$param["Slope",1]*x
  if(scale == "SD")  y <- try(sqrt(y))
  if(scale == "VAR") plot(object$data$x, object$data$y2, main=main, xlab=xlab, ylab=ylab, col=col, ...)
  if(scale == "SD")  plot(object$data$x, sqrt(object$data$y2), main=main, xlab=xlab, ylab=ylab, col=col, ...)
  lines(x, y)
  if(print_param) legend("topleft", legend = 
                           c(as.expression(bquote(italic(a) == .(round_and_format(object$param["a", 1], sign_digits = digits_param)) ~ "\u00B1" ~
                                                    .(round_and_format(object$param["a", 2], sign_digits = digits_param)) ~~~
                                                    italic(b) == .(round_and_format(object$param["b", 1], sign_digits = digits_param)) ~ "\u00B1" ~
                                                    .(round_and_format(object$param["b", 2], sign_digits = digits_param)))),
                             as.expression(bquote(italic(A) == .(round_and_format(object$param["Intercept", 1], sign_digits = digits_param)) ~ "\u00B1" ~ 
                                                    .(round_and_format(object$param["Intercept", 2], sign_digits = digits_param)) ~~~ 
                                                    italic(B) == .(round_and_format(object$param["Slope", 1], sign_digits = digits_param)) ~ "\u00B1" ~ 
                                                    .(round_and_format(object$param["Slope", 2], sign_digits = digits_param)))),
                             as.expression(bquote(italic(R)^2 == .(round_and_format(100*object$Rsquared, digits_rsquared)) ~ "%"))),
                         box.lty = 0, bg="transparent", xjust=0)
}


signif(10000.51, 2)

#' Simulate responses from \code{\link{rate_gls}} fit
#' 
#' \code{simulate.rate_gls} Simulate responses from \code{\link{rate_gls}} fit
#' 
#' \code{simulate.rate_gls} ... 
#' 
#' @param object fitted object from \code{\link{rate_gls}}
#' @param nsim numder of simulations.
#' @return \code{simulate.rate_gls} returns a list where each slot is an interation. 
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @export


simulate.rate_gls <- function(object, nsim = 10){
  sim_out <- list()
  for(i in 1:nsim){
    sim_out[[i]] <-
      try(
        simulate_rate(tree=object$tree, 
                      startv_x = ifelse(object$model == "predictor_gBM", 1, 0),
                      sigma_x = sqrt(object$param["sigma2_x", 1]),
                      a = object$param["a",1],
                      b = object$param["b",1],
                      x = ifelse(rep(object$model == "recent_evol", length(object$data$x)), 
                                 c(object$data$x), NULL),
                      sigma_y = ifelse(object$model == "recent_evol", sqrt(object$sigma2_y), NULL),
                      model = object$model),
      silent = TRUE)
  }
  return(sim_out)
}

#' Boostrap of the rate gls model fit
#' 
#' \code{boot_rate_gls} Boostrap of the rate gls model fit
#' 
#' \code{boot_rate_gls} Provides a parametric boostrap of the generalized least squares rate model fit 
#' using \code{\link{simulate.rate_gls}}  
#' 
#' @param mod output from \code{\link{rate_gls}}.
#' @param n number of bootsrap samples
#' @param useLFO when calulating the mean vector of the traits in the "recent_evol" analysis, should the focal species
#' be left out when calculating the corresponding species' mean. The correct way is to use TRUE, but in practice it has little effect and FALSE 
#' will speed up the model fit (particularely useful when bootstrapping) 

#' @return \code{boot_rate_gls} returns a table with the original estimates and SE from the GLS fit in the two first columns
#' followed by the bootstrap estiamte of the SE and the 2.5\%, 50\% and 97.5\% quantiles of the boostrap distribution. 
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @export

boot_rate_gls <- function(object, n = 10, useLFO = TRUE){
  boot_distribution <- matrix(NA, ncol = 6, nrow = n)
  colnames(boot_distribution) <- c("Intercept", "Slope", "a", "b", "sigma2_x", "Rsquared")
  sim_out <- simulate(object, nsim = n)
  for(i in 1:n){
     mod <- rate_gls(x=sim_out[[i]][,"x"], y=sim_out[[i]][,"y"], species = sim_out[[i]][,"species"], tree=object$tree, 
                     startv = list(object$param["a",1], object$param["b",1]), model = object$model, silent = TRUE,
                     useLFO = useLFO)
    boot_distribution[i,] <- c(mod$param[,1], mod$Rsquared)
    print(paste("Bootstrap iteration", i))
  }
  return(cbind(rbind(object$param, Rsquared = c(object$Rsquared, NA)),
               boot_SE = apply(boot_distribution, 2, function(x) sqrt(var(x))), 
               t(apply(boot_distribution, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975))))))
}

