#' Generalized least squares rate model
#' 
#' \code{rate_gls} Generalized least squares rate model
#' 
#' \code{rate_gls} Fits a generalized least squares model to estimate parameters of the evolutoinary model of two traits x and y, 
#' where the rate of y depends on the value of x. Three models are implmented. In the two first, "predictor_BM" and
#' "predictor_geometricBM", the evolution of y follows a Brownian motion with variance linear in x, while the evolution of x either
#' follows a Brownian motion or a geometric brownian motion, respectively. In the third model, the residuals of the macroveoluitonary
#' predictions of y have variance linear in x. 
#'
#' @param x trait values must be equal to the length of y and tips on the tree. Note that x is mean centred in the analysis
#' @param y trait values
#' @param species names of the species, must be equal in length and in the same order as x and y
#' @param tree object of class \code{\link{phylo}}, needs to be ultrametric and with total length of unit,
#' tips must have the same names as in \code{species}
#' @param model either "predictor_BM", "predictor_geometricBM" or "residual_rate"
#' @param startv vector of starting values for the a and b parameters. ----- NB! for "residual
#' rate" model, the starting value of a is not used. -----
#'  
#' @return \code{rate_gls} 
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

# TODO:
# Check if the startvalue of a should be changed for the residual rate model

rate_gls <- function(x, y, species, tree, model = "predictor_BM", startv = list(a = NULL, b = NULL), maxiter = 100, silent = FALSE){
  ### phylogenetic relatedness matrix ###
  # A <- Matrix::Matrix(ape::vcv(tree), sparse = TRUE)
  # Ainv <- Matrix::solve(A)
  # check out repeated-entry sparce matrix class (by Steve Walker???, is this the "dgTMatrix" class)
  A <- ape::vcv(tree)
  A <- A[match(species, colnames(A)),match(species, colnames(A))] #ordering A according to species
  Ainv <- solve(A)
  
  # if(length(unique(diag(A)) != 1)) stop("The tree is not ultrametric with unit depth")
  # if(unique(diag(A)) == 1)         stop("The tree is not ultrametric with unit depth")
  
  ### predictor (x variable)
  x <- as.matrix(x)
  
  # evolutionary rate
  C <- t(chol(A)) #left cholesky factor
  if(model == "predictor_BM")          sigma2 <- var(solve(C)%*%x) 
  if(model == "predictor_geometricBM") sigma2 <- var(solve(C)%*%log(x))
  if(model == "residual_rate")         sigma2 <- var(x)
  sigma2 <- c(sigma2)
  sigma2_SE <- sqrt(2*sigma2^2/(length(x)+2)) # from Lynch and Walsh 1998 eq A1.10c
  
  # phylogenetic weighted mean
  X <- matrix(rep(1, length(x)), ncol = 1) # design matrix
  if(model == "predictor_BM")          mean_x <- solve(t(X)%*%Ainv%*%X)%*%t(X)%*%Ainv%*%x
  if(model == "predictor_geometricBM") mean_x <- exp(solve(t(X)%*%Ainv%*%X)%*%t(X)%*%Ainv%*%log(x))
  if(model == "residual_rate")         mean_x <- mean(x)
  
  
  ### gls model
  # functions for the different models (to compute a, b, and R)
  if(model == "predictor_BM"){
    a_func <- function(Beta, sigma2, mean_x) Beta[1,1]
    b_func <- function(Beta, sigma2, mean_x) 2*Beta[2,1]
    R_func <- function(a, b, sigma2, t_a, V=NULL) 2*a^2*t_a^2 + c(sigma2*(b/2)^2) * t_a*(1 - 2*t_a + 4*t_a^2)
    a_SE_func <- function(Beta_vcov, sigma2) sqrt(Beta_vcov[1,1])
    b_SE_func <- function(Beta_vcov, sigma2) sqrt(4*Beta_vcov[2,2])
  }
  if(model == "predictor_geometricBM"){
    f <- function(sigma2, t_a = 1){
      X <- try((exp(sigma2*t_a)+2*exp(-0.5*sigma2*t_a)-3)/(exp(sigma2*t_a)-1))
      X[is.na(X)] <- 0
      return(X)
    }
    a_func <- function(Beta, sigma2, mean_x) Beta[1,1]-3*Beta[2,1]*f(sigma2)*(exp(sigma2/2)-1)
    b_func <- function(Beta, sigma2, mean_x) 3*sigma2*Beta[2,1]/(2*f(sigma2))
    R_func <- function(a, b, sigma2, t_a, V=NULL){
      2*a^2*t_a^2 + 8*a*b*t_a/sigma2*(exp(0.5*sigma2*t_a)-1)+
        2*b^2/sigma2^2*(exp(sigma2*t_a)-1)*(exp(sigma2*t_a)-1+2*(exp(0.5*sigma2*t_a)-exp(0.5*sigma2))^2+
                                              4/3*(exp(0.5*sigma2)-exp(0.5*sigma2*t_a))*(f(sigma2, t_a)*exp(0.5*sigma2*t_a)-f(sigma2)*exp(0.5*sigma2))-
                                              4/9*f(sigma2, t_a)*exp(0.5*sigma2*t_a)*f(sigma2)*exp(0.5*sigma2)+
                                              2/9*f(sigma2)^2*exp(sigma2)
        )
    }
    a_SE_func <- function(Beta_vcov, sigma2) sqrt(Beta_vcov[1,1] + 9*f(sigma2)^2*(exp(sigma2/2)-1)^2*Beta_vcov[2,2] - 3*f(sigma2)*(exp(sigma2/2-1)*Beta_vcov[1,2]))
    b_SE_func <- function(Beta_vcov, sigma2) sqrt(9*sigma2^2*Beta_vcov[2,2]/(4*f(sigma2)^2))
  }
  if(model == "residual_rate"){
    a_func <- function(Beta, sigma2, mean_x) (1+mean_x^2/sigma2)*Beta[1,1] - Beta[2,1]*mean_x  #can be simplified due to mean centring of x
    b_func <- function(Beta, sigma2, mean_x) Beta[2,1] - Beta[1,1]*mean_x/sigma2               #can be simplified due to mean centring of x
    R_func <- function(a, b, sigma2, t_a, V){
      Vinv <- solve(V)
      dVmin2 <- diag(x = diag(Vinv)^-2)
      return(2* (dVmin2%*%Vinv) * (dVmin2%*%Vinv))
    }
    a_SE_func <- function(Beta_vcov, sigma2) sqrt(Beta_vcov[1,1]) # NB! x must be mean centred for this to be true
    b_SE_func <- function(Beta_vcov, sigma2) sqrt(Beta_vcov[2,2]) # NB! x must be mean centred for this to be true
  }

  # initial values
  RSS <- NA                                       # residual sum of squares
  t_a <- as.vector(A[!lower.tri(A)])              # age of common ancestor        
  a   <- startv$a                                 # parameter of the process
  b   <- startv$b                                 # parameter of the process
  x <- x-c(mean_x)                                # NB! x is mean centred in the analysis
  mean_x_original <- mean_x
  mean_x <- 0
  
  # response variabe and residual covariance matrix
  if(model != "residual_rate") {
    # response variable
    y2 <- y^2
    
    # finding starting values when not specified
    if(is.null(a) | is.null(b)){ 
      coef_lm <- coef(lm(y2~x))
      if(is.null(a)) a <- coef_lm[1]
      if(is.null(b)) b <- coef_lm[2]
    }
    R <- A                                          
    R[!lower.tri(R)] <- R_func(a, b, sigma2, t_a, V)
    R[lower.tri(R)] <- t(R)[lower.tri(R)]
  }
  
  if(model == "residual_rate") {
    mod_Almer <- Almer(y ~ 1 + (1|species), A = list(species = Matrix::Matrix(A, sparse = TRUE)))
    sigma2_y <- VarCorr(mod_Almer)$species[1]
    residvar <- attr(VarCorr(mod_Almer), "sc")^2
    V <- sigma2_y*A + residvar*diag(length(x))
    Vinv <- solve(V)
    y_mean <- as.matrix(apply(cbind(y), 1, function(x) (sum(y)-x)/(length(y)-1)))
    y_predicted <- y_mean + (V - solve(diag(x = diag(Vinv))))%*%Vinv%*%(y-y_mean)
    y2 <- (y-y_predicted)^2
    
    if(is.null(a) | is.null(b)){ # finding starting values
      coef_lm <- coef(lm(y2~x))
      if(is.null(a)) a <- coef_lm[1]
      if(is.null(b)) b <- coef_lm[2]
      E <- c(a)*diag(length(x)) + diag(x = c(b*x))
      V <- sigma2_y*A + E
      R <- R_func(a, b, sigma2, t_a, V)
      
      # Making sure that the correlations in the residuals are below 1 
      while(max(cov2cor(R)[upper.tri(R)])>=1){
        b <- b*0.95
        E <- c(a)*diag(length(x)) + diag(x = c(b*x))
        V <- sigma2_y*A + E
        R <- R_func(a, b, sigma2, t_a, V)
      }
    }
    else{ # when starting values are specified
      E <- c(a)*diag(length(x)) + diag(x = c(b*x))
      V <- sigma2_y*A + E
      R <- R_func(a, b, sigma2, t_a, V)
    }
  }
  
  # inverse of residual covariance matrix
  Rinv <- solve(R)
  
  # gls design matrix  
  X <- as.matrix(cbind(rep(1, length(x)), x))
  
  # updating values using gls
  for(i in 1:maxiter){
    
    if(model == "residual_rate") {
      # updating response variable updating response variable according to new V
      Vinv <- solve(V)
      y_predicted <- y_mean + (V - solve(diag(x = diag(Vinv))))%*%Vinv%*%(y-y_mean)
      y2 <- (y-y_predicted)^2
    }
    
    # gls estimates
    Beta <- solve(t(X)%*%Rinv%*%X)%*%t(X)%*%Rinv%*%y2
    a <- a_func(Beta, sigma2, mean_x)
    b <- b_func(Beta, sigma2, mean_x)
    
    # updating residual covariance matrix using gls estimates
    if(model != "residual_rate") {
      R[!lower.tri(R)] <- R_func(a, b, sigma2, t_a)
      R[lower.tri(R)] <- t(R)[lower.tri(R)]
    }
    
    if(model == "residual_rate") {
      E <- c(a)*diag(length(x)) + diag(x = c(b*x))
      V <- sigma2_y*A + E
      R <- R_func(a, b, sigma2, t_a, V)
      
    }
    
    Rinv <- solve(R)
    
    # residual sum of squares
    RSS[i] <- t(y2-X%*%Beta)%*%Rinv%*%(y2-X%*%Beta)
    if(!silent) print(paste("Generalized residual sum of squares:", RSS[i]))
    if(i>1){
      if(RSS[i]<=RSS[i-1] & RSS[i-1]-RSS[i]<1e-8) break()
    }
  }
  
  Beta_vcov <- solve(t(X)%*%Rinv%*%X)
  Beta_SE <- sqrt(diag(Beta_vcov))
  param <- cbind(rbind(Beta, a, b, sigma2), rbind(cbind(Beta_SE), a_SE_func(Beta_vcov, sigma2), b_SE_func(Beta_vcov, sigma2), sigma2_SE))
  colnames(param) <- c("Estimate", "SE")
  rownames(param) <- c("Intercept", "Slope", "a", "b", "Sigma^2")
  # TODO: include Sigma2_y and phylo_heritability of y in the parameter output? and maybe also mean x. Make a param_y and a param_x output

  TSS <- t(y2)%*%Rinv%*%(y2) # (generalized) total sum of squares
  Rsquared <- 1-(RSS[i]/TSS)
  
  # V_scaled <- R/diag(R)[1]
  # C <- t(chol(V_scaled))
  # C_inv <- solve(C)
  # resid_var <- sum((C_inv%*%(y2-X%*%Beta))^2)/(length(y2)-nrow(Beta))
  
  if(i == maxiter){
    warning("Model reached maximum number of iterations without convergence")
    convergence <- "No convergence"
  }
  else  convergence <- "Convergence"
  
  return(list(param = param, RSS = RSS, Rsquared = Rsquared, R=R, tree = tree, model = model, convergence = convergence))
}
