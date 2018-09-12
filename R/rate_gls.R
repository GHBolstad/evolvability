#' Generalized least squares rate model
#' 
#' \code{rate_gls} Generalized least squares rate model
#' 
#' \code{rate_gls} Fits a generalized least squares model to estimate parameters of the evolutoinary model of two traits x and y, 
#' where the rate of y depends on the value of x. x evolves eiter accoriding to a Brownian motion model of evolution or a geometric 
#' Brownian motion  model of evolution. 
#'
#' @param x trait values must be equal to the length of y and tips on the tree 
#' @param y trait values 
#' @param model either "predictor_BM" or "predictor_geometricBM"
#' @param Beta starting values for the a nad b parameters 
#'  
#'    
#' 
#' @importFrom ape vcv


rate_gls <- function(x, y, tree, model = "predictor_BM", Beta = as.matrix(c(mean(y^2), 0)), maxiter = 100, silent = FALSE){
  ### phylogenetic relatedness matrix ###
  # A <- Matrix::Matrix(ape::vcv(tree), sparse = TRUE)
  # Ainv <- Matrix::solve(A)
  # check out repeated-entry sparce matrix class (by Steve Walker???, is this the "dgTMatrix" class)
  A <- ape::vcv(tree)
  Ainv <- solve(A)
  
  # if(length(unique(diag(A)) != 1)) stop("The tree is not ultrametric with unit depth")
  # if(unique(diag(A)) == 1)         stop("The tree is not ultrametric with unit depth")
  
  ### predictor (x variable)
  # evolutionary rate
  C <- t(chol(A)) #left cholesky factor
  if(model == "predictor_BM")          sigma2 <- var(solve(C)%*%x) 
  if(model == "predictor_geometricBM") sigma2 <- var(solve(C)%*%log(x))
  if(model == "residual_rate")         sigma2 <- var(x)
  # sigma2_SE <- sqrt(2*sigma2^2/(length(x)+2)) # from Lynch and Walsh 1998 eq A1.10c
  
  # phylogenetic weighted mean
  X <- matrix(rep(1, length(x)), ncol = 1) # design matrix
  V <- A*c(sigma2) # residual covariance matrix
  if(model == "predictor_BM")          mean_x <- solve(t(X)%*%Ainv%*%X)%*%t(X)%*%Ainv%*%x
  if(model == "predictor_geometricBM") mean_x <- exp(solve(t(X)%*%Ainv%*%X)%*%t(X)%*%Ainv%*%log(x))
  if(model == "residual_rate")         mean_x <- mean(x)
  
  
  ### gls model
  
  # functions for the different models (to compute a, b, and V)
  if(model == "predictor_BM"){
    a_func <- function(Beta, sigma2) Beta[1,1]
    b_func <- function(Beta, sigma2) 2*Beta[2,1]
    V_func <- function(a, b, sigma2, t_a) 2*a^2*t_a^2 + c(sigma2*(b/2)^2) * t_a*(1 - 2*t_a + 4*t_a^2)
    #a_SE_func <- function(Beta_vcov, sigma2) sqrt(Beta_vcov[1,1])
    #b_SE_func <- function(Beta_vcov, sigma2) sqrt(4*Beta_vcov[2,2])
  }
  if(model == "predictor_geometricBM"){
    f <- function(sigma2, t_a = 1){
      X <- try((exp(sigma2*t_a)+2*exp(-0.5*sigma2*t_a)-3)/(exp(sigma2*t_a)-1))
      X[is.na(X)] <- 0
      return(X)
    }
    a_func <- function(Beta, sigma2) Beta[1,1]-3*Beta[2,1]*f(sigma2)*(exp(sigma2/2)-1)
    b_func <- function(Beta, sigma2) 3*sigma2*Beta[2,1]/(2*f(sigma2))
    V_func <- function(a, b, sigma2, t_a){
      2*a^2*t_a^2 + 8*a*b*t_a/sigma2*(exp(0.5*sigma2*t_a)-1)+
        2*b^2/sigma2^2*(exp(sigma2*t_a)-1)*(exp(sigma2*t_a)-1+2*(exp(0.5*sigma2*t_a)-exp(0.5*sigma2))^2+
                                              4/3*(exp(0.5*sigma2)-exp(0.5*sigma2*t_a))*(f(sigma2, t_a)*exp(0.5*sigma2*t_a)-f(sigma2)*exp(0.5*sigma2))-
                                              4/9*f(sigma2, t_a)*exp(0.5*sigma2*t_a)*f(sigma2)*exp(0.5*sigma2)+
                                              2/9*f(sigma2)^2*exp(sigma2)
        )
    }
    #a_SE_func <- function(Beta_vcov, sigma2) sqrt(Beta_vcov[1,1] + 9*f(sigma2)^2*(exp(sigma2/2)-1)^2*Beta_vcov[2,2] - 3*f(sigma2)*(exp(sigma2/2-1)*Beta_vcov[1,2]))
    #b_SE_func <- function(Beta_vcov, sigma2) sqrt(9*sigma2^2*Beta_vcov[2,2]/(4*f(sigma2)^2))
  }
  
  # initial values
  RSS <- NA                                       # residual sum of squares
  t_a <- as.vector(A[!lower.tri(A)])              # age of common ancestor        
  a   <- as.vector(a_func(Beta, sigma2))          # parameter of the process
  b   <- as.vector(b_func(Beta, sigma2))          # parameter of the process
  V <- A                                          # residual covariance matrix
  V[!lower.tri(V)] <- V_func(a, b, sigma2, t_a)
  V[lower.tri(V)] <- t(V)[lower.tri(V)]
  Vinv <- solve(V)
  
  # response variable
  if(model == "residual_rate") {
    y_mean <- as.matrix(apply(cbind(y), 1, function(x) (sum(y)-x)/(length(y)-1)))
    y_predicted <- y_mean +(V - solve(diag(x = diag(Vinv))))%*%Vinv%*%(y-y_mean)
    y2 <- (y-y_predicted)^2
  }
  else y2 <- y^2
  
  # design matrix  
  X <- as.matrix(cbind(rep(1, length(x)), x))  ## x - mean(x) have changed so that x is not mean centred
  
  # updating values
  for(i in 1:maxiter){
    Beta <- solve(t(X)%*%Vinv%*%X)%*%t(X)%*%Vinv%*%y2
    a <- a_func(Beta, sigma2)
    b <- b_func(Beta, sigma2)
    V[!lower.tri(V)] <- V_func(a, b, sigma2, t_a)
    V[lower.tri(V)] <- t(V)[lower.tri(V)]
    Vinv <- solve(V)
    
    # residual sum of squares
    RSS[i] <- t(y2-X%*%Beta)%*%Vinv%*%(y2-X%*%Beta)
    if(!silent) print(paste("Residual sum of squares:", RSS[i]))
    if(i>1){
      if(RSS[i]<=RSS[i-1] & RSS[i-1]-RSS[i]<1e-8) break()
    }
  }
  
  # Beta_vcov <- solve(t(X)%*%Vinv%*%X)
  # Beta_SE <- sqrt(diag(Beta_vcov))
  # param <- cbind(rbind(Beta, a, b, sigma2), rbind(cbind(Beta_SE), a_SE_func(Beta_vcov, sigma2), b_SE_func(Beta_vcov, sigma2), sigma2_SE))
  # colnames(param) <- c("Estimate", "SE")
  # rownames(param) <- c("Intercept", "Slope", "a", "b", "Sigma2")
  
  param <- cbind(rbind(Beta, a, b, mean_x, sigma2))
  colnames(param) <- c("Estimate")
  rownames(param) <- c("Intercept", "Slope", "a", "b", "mean_x", "Sigma^2")
  
  TSS <- t(y2)%*%Vinv%*%(y2) # (generalized) total sum of squares
  Rsquared <- 1-(RSS[i]/TSS)
  
  # V_scaled <- V/diag(V)[1]
  # C <- t(chol(V_scaled))
  # C_inv <- solve(C)
  # resid_var <- sum((C_inv%*%(y2-X%*%Beta))^2)/(length(y2)-nrow(Beta))
  
  if(i == maxiter) print(paste("Model reached maximum number of iterations without convergence"))
  return(list(param = param, RSS = RSS, Rsquared = Rsquared, V=V, tree = tree, model = model))
}
