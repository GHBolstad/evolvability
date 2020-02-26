#' Generalized least squares rate model
#' 
#' \code{rate_gls} fits a generalized least squares model to estimate parameters of 
#' the evolutoinary model of two traits x and y, where the evolutionary rate of y 
#' depends on the value of x. Three models are implmented. In the two first, 
#' "predictor_BM" and "predictor_gBM", the evolution of y follows a Brownian motion 
#' with variance linear in x, while the evolution of x either follows a Brownian 
#' motion or a geometric brownian motion, respectively. In the third model, the 
#' residuals of the macroveoluitonary predictions of y have variance linear in x. 
#'
#' @param x explanatory variable must be equal to the length of y and tips on the 
#' tree. Note that the algorithm mean centres x in the "predictor_BM" and "recent_evol" 
#' analyses, while it is mean stadardized (divided by the mean) in the "predictor_gBM".
#' @param y trait values of response variable. Note that the algorithm mean centres y.
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
#' @param useLFO acronym for use "Leave Focal Out" when calulating the mean vector of 
#' the traits. In the "recent_evol" analysis, should the focal species be left out 
#' when calculating the corresponding species' mean. The correct way is to use TRUE, 
#' but in practice it has little effect and FALSE will speed up the model fit 
#' (particularely useful when bootstrapping) 
#' 
#' @details The generalized least squares (GLS) model fit a regression where 
#' explanatory variable is x (mean centred) and the response variable is vector of 
#' squared y values for the "predictor_BM" and "predictor_gBM" models and squared 
#' deviation from the evolutionary predictions (squared residuals) in the 
#' "recent_evol" model. From the intercept and slope of the GLS fit the following 
#' evolutionary paramters (a and b) are estimated. The evolutionary models are
#' described in a vignette and in Hansen et al. (in prep)
#' 
#' \deqn{dy = \sqrt{a + bx}dW_{i}}
#' 
#' @return \code{rate_gls} returns a list with elements
#' \itemize{
#'  \item{\code{model}: name of fitted model ("predictor_BM", "predictor_gBM" or "recent_evolution").}
#'  \item{\code{param}: parameter estimates and standard errors, where "A" and 
#'  "B" are the intercept and slope of the GLS regression, "a" and "b" are parameters 
#'  of the evolutionary models (see vignette), and "s2" is the BM-rate parameter of x for the "predictor_BM" model,
#'  the BM-rate parameter for log x for the "predictor_gBM" model, and the variance of x for the 
#'  "recent_evolution" model.}
#'  \item{\code{R2}: the generalized R squared of the GLS model fit}
#'  \item{GLS_objective_score}{The score of the GLS objective function}
#'  \item{R}{Residual variance matrix}
#'  \item{tree}{The phylogenetic tree}
#'  \item{data}{The data used for the GLS regresssion}
#'  \item{convergence}{Whether the algoritm converged or not}
#'  \item{call}
#'  }
#'  
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @importFrom ape vcv
#' @importFrom lme4 VarCorr
#' @importFrom Matrix Matrix solve Diagonal rowSums
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

  #### storing original y and x values ####
  y_original <- y
  x_original <- x
  
  #### y-variable ####
  # mean (weighted assuming BM process and the phylogeny)
  X <- matrix(rep(1, length(y)), ncol = 1) # design matrix
  mod <- GLS(y, X, R = A)$coef
  mean_y <- mod[1]
  Vy <- mod[2]^2
  
  # mean centring
  y <- y - c(mean_y) # note that for the residual rate model, this does not have any effect as the squared deviation from the predictons are used
  #Vy <- 0 # NB!!!!! must be removed when mean scaling !!!!!!
  #L <- t(chol(A))        #!!!!!!!!!CHOLESKY TEST!!!!!!!!
  #y <- forwardsolve(L, y)#!!!!!!!!!CHOLESKY TEST!!!!!!!!
  
  #### x-variable ####
  if(model == "predictor_BM"){
    mod <- GLS(x, X, A)
    mean_x <- mod$coef[1]
    Vx <- mod$coef[2]^2
    s2 <- mod$sigma2
    #s2 <- 0.25^2 # Setting it to true value!
    x <- x - c(mean_x)
    #Vx <- 0 # NB!!!!! must be removed when mean scaling !!!!!!
    AoA <- A*A
    #I <- diag(nrow(A))
    x <- c(x - (1/2)*AoA%*%solve(A, x)) # solve(A, x) equals solve(A)%*%x, but is numerically more stable (and faster)
    #Ainv <- Matrix::solve(A)
    #x <- as.matrix((I-(1/2)*AoA%*%Ainv)%*%x)
    #x <- forwardsolve(L, x) #!!!!!!!!!CHOLESKY TEST!!!!!!!!
  }
  if(model == "predictor_gBM"){
    mod <- GLS(log(x), X, A)
    mean_x <- exp(mod$coef[1])
    Vx <- mod$coef[2]^2
    s2 <- mod$sigma2
    #s2 <- 2 # Setting it to true value!
    if(s2 < 1) warning("sigma(x)^2 < 1, the estimates are innaccurate, probably better to use model = 'predictor_BM'")
    x <- x/c(mean_x)
    #x <- x-1
    #Vx <- 0 # NB!!!!! must be removed when mean scaling !!!!!!
    e_32s2A   <- exp(3/2*s2*A)-1
    e_s2A <- exp(s2*A)-1
    x <- c((2/s2)*(x-(2/3)*exp(-s2/2)*e_32s2A%*%solve(e_s2A, x)))  
    #x <- x-mean(x) ########REMOVE THIS!!!!!
    #I <- diag(nrow(A))
    #inv_e_s2A <- solve(exp(s2*A)-1)
    #x <- as.matrix((2/s2)*(I-(2/3)*exp(-s2/2)*e_32s2A%*%inv_e_s2A)%*%x) 
  }
  if(model == "recent_evol"){
    mean_x <- mean(x)
    s2 <- var(x)
    Vx <- s2/length(x)
    x <- x - c(mean_x) 
  }
  s2 <- c(s2)
  s2_SE <- sqrt(2*s2^2/(length(x)+2)) # from Lynch and Walsh 1998 eq A1.10c 

  #A <- diag(length(species)) #!!!!!!!!!CHOLESKY TEST!!!!!!!!
  
  
  
#### Internal functions ####
  if(model == "predictor_BM"){
    a_func <- function(Beta, a_bias) Beta[1,1] + Vy
    b_func <- function(Beta) Beta[2,1]
    Q <- (A*A*A) - (1/4)*AoA%*%solve(A, AoA) # solve(A, AoA) is equivalent to solve(A)%*%AoA
    #Q <- (A*A*A) - (1/4)*AoA%*%Ainv%*%AoA #by pulling this outside the function, it does not need to be calculated each iteration
    R_func <- function(a, b) 4*a*Vy*A + 2*(a^2 + b^2*Vx)*AoA + b^2*s2*Q
    a_SE_func <- function(Beta_vcov) sqrt(Beta_vcov[1,1]) #CAN WE GET THE UNCERTAINTY IN Vy????? approximation using Lynch and Walsh formula?
    b_SE_func <- function(Beta_vcov) sqrt(Beta_vcov[2,2])
  }
  if(model == "predictor_gBM"){
    # f <- function(s2, t_a = 1){
    #   X <- try((exp(s2*t_a)+2*exp(-0.5*s2*t_a)-3)/(exp(s2*t_a)-1))
    #   X[t_a == 0] <- 0
    #   return(X)
    # }
    # 
    # b_func <- function(Beta, s2) 3*s2*Beta[2,1]/(2*f(s2))
    # #a_func <- function(Beta, s2, a_bias) Beta[1,1]-3*Beta[2,1]*f(s2)*(exp(s2/2)-1)
    # a_func <- function(Beta, s2, a_bias) Beta[1,1]-2*b_func(Beta, s2)*(exp(s2/2)-1)/s2
    # R_func <- function(a, b, s2, t_a){
    #   2*a^2*t_a^2 + (8*a*b*t_a/s2)*(exp(0.5*s2*t_a)-1)+
    #     (2*b^2/s2^2)*(exp(s2*t_a)-1)*(exp(s2*t_a)-1+2*(exp(0.5*s2*t_a)-exp(0.5*s2))^2+
    #                                               4/3*(exp(0.5*s2)-exp(0.5*s2*t_a))*(f(s2, t_a)*exp(0.5*s2*t_a)-f(s2)*exp(0.5*s2))-
    #                                               4/9*f(s2, t_a)*exp(0.5*s2*t_a)*f(s2)*exp(0.5*s2)+
    #                                               2/9*f(s2)^2*exp(s2))
    #     
    # }
    # a_SE_func <- function(Beta_vcov, s2) sqrt(Beta_vcov[1,1] + 9*f(s2)^2*(exp(s2/2)-1)^2*Beta_vcov[2,2] - 3*f(s2)*(exp(s2/2-1)*Beta_vcov[1,2]))
    # b_SE_func <- function(Beta_vcov, s2) sqrt(9*s2^2*Beta_vcov[2,2]/(4*f(s2)^2))
    
    a_func <- function(Beta, a_bias) Beta[1,1] + Vy - 2*Beta[2,1]*exp(Vx/2)*(exp(s2/2) - 1)/s2
    #a_func <- function(Beta, a_bias) Beta[1,1] + Vy -  2*0.1      *exp(Vx/2)*(exp(s2/2) - 1)/s2
    b_func <- function(Beta) Beta[2,1]
    # some matrix calculations to avoid repeating them in the loop ...
    AoA <- A*A
    e_hlfVx  <- exp(Vx/2)
    e_2Vx    <- exp(2*Vx)
    e_hlfs2A <- exp(s2/2*A)-1
    e_2s2A   <- exp(2*s2*A)-1
    Q1       <- 8*e_hlfVx/s2*A*e_hlfs2A 
    Q2       <- 2*e_2Vx/s2^2*(8/3*(e_2s2A-e_hlfs2A)-e_2s2A-8/9*e_32s2A%*%solve(e_s2A, e_32s2A, tol = 1e-99))
    # ... end matrix calculations
    R_func <- function(a, b) 4*a*Vy*A + 8*b*Vy*e_hlfVx/s2*e_hlfs2A + 2*a^2*AoA + a*b*Q1 + b^2 * Q2
    a_SE_func <- function(Beta_vcov) sqrt(Beta_vcov[1,1]) # This is only an apprixmation, check if it good with parametric simulations
    b_SE_func <- function(Beta_vcov) sqrt(Beta_vcov[2,2])
  }
  if(model == "recent_evol"){
    a_func <- function(Beta, a_bias) Beta[1,1]-a_bias  
    b_func <- function(Beta) Beta[2,1]
    R_func <- function(V, V_micro){
      Vinv <- solve(V)
      dVinv <- diag(x = diag(Vinv))
      Q <- solve(dVinv%*%V%*%dVinv)
      I <- diag(nrow(A))
      a_bias <- mean(diag(Q + (I-2*solve(dVinv)%*%Vinv)%*%V_micro))
      return(list(R = 2*Q*Q, a_bias = a_bias))
    }
    a_SE_func <- function(Beta_vcov) sqrt(Beta_vcov[1,1]) # !!!! This is wrong as it assumes that a_bias is measured without uncertainty !!!! (check if it is a good approximation)
    b_SE_func <- function(Beta_vcov) sqrt(Beta_vcov[2,2])
  }
  
  
  
  #~~~~~~~~~~~~~~~~~~~~~#
  #### The GLS model ####
  #~~~~~~~~~~~~~~~~~~~~~#
  
  #### Response variabe and initial values ####
  obj <- NULL                                     # objective function scores
  a   <- startv$a                                 # parameter of the process
  b   <- startv$b                                 # parameter of the process
  
  if(model == "predictor_BM" | model == "predictor_gBM"){
    # age of common ancestor
    t_a <- as.vector(A[which(!lower.tri(A))])     #  not needed anymore...        
    
    # response variable
    y2 <- y^2

    # finding starting values when not specified
    if(is.null(a) | is.null(b)){ 
      coef_lm <- coef(lm(y2~x))
      if(is.null(a)) a <- coef_lm[1]
      if(is.null(b)) b <- 2*coef_lm[2]
    }
    R <- matrix(nrow = nrow(A), ncol = ncol(A))                                    
  }
  if(model == "recent_evol") {
    mod_Almer <- Almer(y ~ 1 + (1|species), A = list(species = Matrix::Matrix(A, sparse = TRUE)))
    #mod_Almer <- Almer(y ~ 1 + (1|species), A = list(species = A))
    sigma2_y <- lme4::VarCorr(mod_Almer)$species[1]
    sigma2_y_SE <- sqrt(2*sigma2_y^2/(length(y)+2)) # from Lynch and Walsh 1998 eq A1.10c 
    residvar <- attr(lme4::VarCorr(mod_Almer), "sc")^2
    if(is.null(a)) a <- residvar 
    if(is.null(b)) b <- 0
  }
  
  # GLS design matrix  
  X <- as.matrix(cbind(rep(1, length(x)), x))
  a_bias <- NULL # this is only needed for the "recent_evol" model
  
  #### Iterative GLS ####
  for(i in 1:maxiter){
    if(model == "predictor_BM" | model == "predictor_gBM"){
      R <- R_func(a, b) # R_func(1, 0.1) #
      #R[!lower.tri(R)] <- R_func(a, b, s2, t_a)
      #R[lower.tri(R)] <- t(R)[lower.tri(R)]
    }
    
    if(model == "recent_evol") {
      diag_V_micro <- a + b*x
      #diag_V_micro[diag_V_micro < 0] <- 0 # Negative residual variances are replaced by zero
      V_micro <- diag(diag_V_micro)
      V <- sigma2_y*A + V_micro
      R_list <- R_func(V, V_micro)
      R <- R_list$R
      a_bias <- R_list$a_bias
      y_predicted <- macro_pred(y=y, V=V, useLFO=useLFO) # maybe change this to using equation 14 instead
      y2 <- (y-y_predicted)^2
    }
    
    # GLS estimates
    mod <- GLS(y=y2, X, R, coef_only = TRUE) 
    Beta <- cbind(mod$coef)
    a <- a_func(Beta, a_bias)
    b <- b_func(Beta)
    minimum_a <- c(-min(b*x))
    if(model == "recent_evol" & a < minimum_a) a <- minimum_a # forces the microevolutionary variance to be zero or higher
    
    # Objective function
    obj[i] <- mod$GSSE
    if(!silent) print(paste("Generalized sum of squares:", obj[i]))
    if(i>1){
      if(obj[i]<=obj[i-1] & obj[i-1]-obj[i]<1e-7) break()
    }
  }
  
  
  #### CAN BE DELETED!!! TEST's THAT THE GLS GIVE WHAT IT SHOULD####
  # Rtest <- R
  # GLS(y=y2, X, Rtest, coef_only = TRUE)
  # param <- solve(t(X)%*%solve(Rtest)%*%X)%*%t(X)%*%solve(Rtest)%*%y2
  # param
  # t(y2-X%*%param)%*%solve(Rtest)%*%(y2-X%*%param)
  ##############
  
  mod <- GLS(y2, X, R) # rerunning the model with all the output
  Beta_vcov <- mod$coef_vcov
  #Beta_SE <- sqrt(diag(Beta_vcov))
  param <- cbind(rbind(a, b, s2), rbind(a_SE_func(Beta_vcov), b_SE_func(Beta_vcov), s2_SE))
  colnames(param) <- c("Estimate", "SE")
  rownames(param) <- c("a", "b", "sigma(x)^2") 
  # REMOVE A and B altoghether (maybe A is neede for plottin, or recalculate it in the plotting function)
  # change sigma(x)^2 to var(x) ??? should be easier to understand 

  if(model == "recent_evol"){
    param <- rbind(param, cbind(sigma2_y, sigma2_y_SE) )
    rownames(param)[4] <- "sigma(y)^2"
  }
  
  # R squared
  Rsquared <- mod$R2

  
  if(i == maxiter){
    warning("Model reached maximum number of iterations without convergence")
    convergence <- "No convergence"
  }
  else  convergence <- "Convergence"
  
  report <-  list(model = model, param = param, Rsquared = Rsquared, GLS_objective_scores = obj,
                  R=R, Beta = Beta, Beta_vcov = Beta_vcov, tree = tree, data = list(y2 = y2, x = x, y = y), convergence = convergence,
                  additional_param = c(mean_y = mean_y, Vy = Vy, mean_x = mean_x, Vx = Vx))
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
#' @param scale the scale of the y-axis, either the variance scale ("VAR"), that is y^2, or the standard deviation scale ("SD"), that is abs(y).
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
  y <- object$Beta[1] + object$Beta[2]*x
  if(scale == "SD")  y <- try(sqrt(y))
  if(scale == "VAR") plot(object$data$x, object$data$y2, main=main, xlab=xlab, ylab=ylab, col=col, ...)
  if(scale == "SD")  plot(object$data$x, sqrt(object$data$y2), main=main, xlab=xlab, ylab=ylab, col=col)#, ...)
  lines(x, y)

  if(print_param) legend("topleft", legend = 
                           c(as.expression(bquote(italic(a) == .(round_and_format(object$param["a", 1], sign_digits = digits_param)) ~ "\u00B1" ~
                                                    .(round_and_format(object$param["a", 2], sign_digits = digits_param)) ~~~
                                                    italic(b) == .(round_and_format(object$param["b", 1], sign_digits = digits_param)) ~ "\u00B1" ~
                                                    .(round_and_format(object$param["b", 2], sign_digits = digits_param)))),
                             # as.expression(bquote(italic(A) == .(round_and_format(object$param["A", 1], sign_digits = digits_param)) ~ "\u00B1" ~ 
                             #                        .(round_and_format(object$param["A", 2], sign_digits = digits_param)) ~~~ 
                             #                        italic(B) == .(round_and_format(object$param["B", 1], sign_digits = digits_param)) ~ "\u00B1" ~ 
                             #                        .(round_and_format(object$param["B", 2], sign_digits = digits_param)))),
                             as.expression(bquote(italic(R)^2 == .(round_and_format(100*object$Rsquared, digits_rsquared)) ~ "%"))),

                         box.lty = 0, bg="transparent", xjust=0)
}



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
                      sigma_x = sqrt(object$param["sigma(x)^2", 1]),
                      a = object$param["a",1],
                      b = object$param["b",1],
                      x = ifelse(rep(object$model == "recent_evol", length(object$data$x)), c(object$data$x), NULL),
                      sigma_y = ifelse(object$model == "recent_evol", sqrt(object$param["sigma(y)^2", 1]), NULL),
                      model = object$model),
        silent = TRUE)
  }
  return(sim_out)
}

#' Boostrap of the rate gls model fit
#' 
#' \code{boot.rate_gls} Boostrap of the rate gls model fit
#' 
#' \code{boot_rate_gls} Provides a parametric boostrap of the generalized least squares rate model fit 
#' using \code{\link{simulate.rate_gls}}  
#' 
#' @param mod output from \code{\link{rate_gls}}.
#' @param n number of bootsrap samples
#' @param useLFO when calulating the mean vector of the traits in the "recent_evol" analysis, should the focal species
#' be left out when calculating the corresponding species' mean. The correct way is to use TRUE, but in practice it has little effect and FALSE 
#' will speed up the model fit (particularely useful when bootstrapping) 
#' @param silent 
#' 
#' @return \code{boot.rate_gls} returns list where the first slot is a table with the original estimates and SE from the GLS fit in the two first columns
#' followed by the bootstrap estiamte of the SE and the 2.5\%, 50\% and 97.5\% quantiles of the boostrap distribution. The second slot is the complete distribution.
#' 
#' @author Geir H. Bolstad
#' 
#' @examples
#' 
#' @export

boot.rate_gls <- function(object, n = 10, useLFO = TRUE, silent = FALSE){
  if(object$model == "recent_evol"){
    boot_distribution <- matrix(NA, ncol = 5, nrow = n)  
    colnames(boot_distribution) <- c("a", "b", "sigma2_x", "sigma2_y", "Rsquared")
  } else{
    boot_distribution <- matrix(NA, ncol = 4, nrow = n)
    colnames(boot_distribution) <- c("a", "b", "s2", "Rsquared")
  }
  for(i in 1:n){
    sim_out <- simulate.rate_gls(object, nsim = 1)
    mod <- rate_gls(x=sim_out[[1]][,"x"], y=sim_out[[1]][,"y"], species = sim_out[[1]][,"species"], tree=object$tree, 
                    startv = list(object$param["a",1], object$param["b",1]), model = object$model, silent = TRUE,
                    useLFO = useLFO)
    boot_distribution[i,] <- c(mod$param[,1], mod$Rsquared)
    if(mod$convergence != "Convergence") boot_distribution[i,] <- NA
    if(!silent) print(paste("Bootstrap iteration", i))
  }
  return(list(summary = cbind(rbind(object$param, Rsquared = c(object$Rsquared, NA)),
                        boot_mean = apply(boot_distribution, 2, mean, na.rm = TRUE),
                        boot_SE = apply(boot_distribution, 2, function(x) sqrt(var(x, na.rm = TRUE))), 
                        t(apply(boot_distribution, 2, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)))),
              boot_distribution = boot_distribution))
}

