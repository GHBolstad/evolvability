## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
tree <- geiger::sim.bdtree(b = 1, d = 0, n = 100, t = 4)
ape::is.ultrametric(tree)

## ------------------------------------------------------------------------
tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1]

## ------------------------------------------------------------------------
sim_data <- simulate_rate(tree, startv_x=0, sigma_x=1, a=2, b=1, model = "predictor_BM")

## ------------------------------------------------------------------------
gls_mod <- rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "predictor_BM", maxiter = 100, silent = TRUE)
gls_mod$param

## ------------------------------------------------------------------------
boot_rate_gls(gls_mod, n = 5) # when doing a proper boostrap n should be larger n>999

## ------------------------------------------------------------------------
plot(gls_mod, scale = "VAR")

## ------------------------------------------------------------------------
plot(gls_mod) # with the default scale == "SD"

## ------------------------------------------------------------------------
sim_data <- simulate_rate(tree, startv_x=0, sigma_x=1, a=2, b=1, sigma_y = 2, model = "recent_evol")

## ------------------------------------------------------------------------
gls_mod <- rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "recent_evol", maxiter = 1000, silent = TRUE, useLFO = FALSE)
gls_mod$param

## ------------------------------------------------------------------------
boot_rate_gls(gls_mod, n = 5, useLFO = FALSE) 
# when doing a proper boostrap n should be larger n>999


## ------------------------------------------------------------------------
plot(gls_mod, scale = "VAR")

## ---- eval = FALSE-------------------------------------------------------
#  
#  
#  #### recent_evol ####
#  # Testing the estimation of b
#  b_true <- c(1:100)/50
#  b_est <- NA
#  
#  for(i in 1:length(b_true)){
#    sim_data <- try(simulate_rate(tree, startv_x=0, sigma_x=3, a=10, b=b_true[i], sigma_y = 2,
#                                  model = "recent_evol"), TRUE)
#    if(is.null(nrow(sim_data))){
#      b_est[i] <- NA
#    }else{
#      gls_mod <- try(rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "recent_evol",
#                          maxiter = 1000, silent = TRUE, use_leave_focal_out_for_y_mean = FALSE), TRUE)
#      if(length(gls_mod) == 1){
#        b_est[i] <- NA
#      }else{
#        b_est[i] <- gls_mod$param["b",1]
#      }
#    }
#  }
#  plot(b_true, b_est)
#  abline(0,1)
#  #lm(b_est~b_true) # this is good
#  
#  # Testing the estimation of a
#  a_true <- c(1:100)/50
#  a_est <- NA
#  for(i in 1:length(a_true)){
#    sim_data <- try(simulate_rate(tree, startv_x=0, sigma_x=3, a=a_true[i], b=0.1, sigma_y = 2,
#                                  model = "recent_evol"), TRUE)
#    if(is.null(nrow(sim_data))){
#      a_est[i] <- NA
#    }else{
#      gls_mod <- try(rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "recent_evol",
#                          maxiter = 1000, silent = TRUE, use_leave_focal_out_for_y_mean = FALSE), TRUE)
#      if(length(gls_mod) == 1){
#        a_est[i] <- NA
#      }else{
#        a_est[i] <- gls_mod$param["a",1]
#      }
#    }
#  }
#  plot(a_true, a_est)
#  abline(0,1)
#  #lm(a_est~a_true) # this is good
#  

## ---- eval = FALSE-------------------------------------------------------
#  # testing the R2 values
#  
#  # star phylogeny:
#  tree <- phytools::starTree(1:100, rep(1, 100))
#  
#  # predictor_BM model
#  mod <- simulate_rate(tree, startv_x=0, sigma_x=1, a=2, b=1, model = "predictor_BM")
#  gls_mod <- rate_gls(x=mod$x, y=mod$y, species=tree$tip.label, tree, model = "predictor_BM", maxiter = 100, silent = FALSE)
#  gls_mod$Rsquared
#  summary(lm(gls_mod$data$y2~gls_mod$data$x))$r.squared
#  
#  plot(gls_mod)
#  
#  X <- simulate(gls_mod, nsim = 10)
#  
#  
#  # predictor_geometricBM model
#  mod <- rate_sim(tree, startv_x=1, sigma_x=1, a=2, b=1, model = "predictor_geometricBM")
#  gls_mod <- rate_gls(x=mod$x, y=mod$y, species=tree$tip.label, tree, model = "predictor_geometricBM", maxiter = 100, silent = FALSE)
#  gls_mod$Rsquared
#  summary(lm(gls_mod$data$y2~gls_mod$data$x))$r.squared
#  
#  # residual_rate model
#  # This one should not be exactly equal
#  mod <- rate_sim(tree, startv_x=0, sigma_x=1, a=2, b=1, model = "residual_rate")
#  gls_mod <- rate_gls(x=mod$x, y=mod$y, species=tree$tip.label, tree, model = "residual_rate", maxiter = 100, silent = FALSE)
#  gls_mod$Rsquared
#  summary(lm(gls_mod$data$y2~gls_mod$data$x))$r.squared
#  
#  # R2 in ideal situation:
#  n <- 100000  #total sample size
#  x <- rep(c(0,1), each = n/2)
#  y <- c(rnorm(n/2, 0, 1), rnorm(n/2, 0, sqrt(2)))
#  summary(lm(y^2 ~ x))$r.squared
#  
#  plot(x, y^2)
#  
#  tree <- phytools::starTree(1:n, rep(1, n))
#  gls_mod <- rate_gls(x=x, y=y, species=tree$tip.label, tree, model = "predictor_BM", maxiter = 100, silent = FALSE)
#  gls_mod$Rsquared
#  summary(lm(gls_mod$data$y2~gls_mod$data$x))$r.squared
#  # NB this is half of what Thomas derived analytically
#  
#  # analytic derivation:
#  # within group var for y^2 is 2*var(y)^2 (derivated from results in Bohrnstedt and Goldberger 1969 + fourth moment of normal = 3*var[y]^2)
#  # avarage within group variation is
#  mean(c(2*1, 2*(2^2))) # = 5
#  # mean of the two grups are 1 and sqrt(2)
#  # the between group mean square is
#  ((2-1)/2)^2 # = 0.25
#  # R2 is
#  1 - 5/(5+0.25)
#  1/21
#  

## ---- eval = FALSE-------------------------------------------------------
#  
#  # testing:
#  tree <- geiger::sim.bdtree(b = 1, d = 0, n = 100, t = 4)
#  tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1]
#  
#  mod <- rate_sim(tree, startv_x=3, sigma_x=1, a= 0.001, b=1, model = "predictor_geometricBM")
#  
#  #mod$y <- mod$y + rnorm(length(mod$y))
#  
#  
#  gls_mod <- rate_gls(x=mod$x, y=mod$y, species=tree$tip.label, tree, model = "predictor_geometricBM", maxiter = 100, silent = FALSE)
#  
#  
#  gls_mod$param["a",1]
#  
#  rate_sim()
#  
#  
#  A <-ape::vcv(tree)
#  C <- t(chol(A))
#  R <- diag(ncol(A))*1
#  
#  sigma_x <- var(solve(C)%*%mod$x)
#  
#  Ainv <- solve(A)
#  Rinv <- solve(R)
#  
#  x1 <- solve(C)%*%(mod$x-mean(mod$x))
#  x2 <- solve(Rinv+solve(c(sigma_x)*A))%*%Rinv%*%(mod$x-mean(mod$x))
#  x3 <- (A-solve(diag(x = diag(Ainv))))%*%Ainv%*%(mod$x-mean(mod$x))
#  x4 <- A%*%solve(A+R)%*%(mod$x-mean(mod$x))
#  
#  plot(x2, x4)
#  plot(mod$x, x4)
#  
#  plot(solve(C)%*%mod$x, solve(Rinv+solve(c(sigma_x)*A))%*%Rinv%*%(mod$x-mean(mod$x)))
#  
#  plot(solve(C)%*%mod$x, (A-solve(diag(x = diag(Ainv))))%*%Ainv%*%(mod$x-mean(mod$x)))
#  
#  
#  #### testing almer
#  tree <- geiger::sim.bdtree(b = 1, d = 0, n = 100, t = 4)
#  tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1]
#  mod <- rate_sim(tree, startv_x=3, sigma_x=1, a= 0.001, b=1)
#  A <- Matrix::Matrix(ape::vcv(tree), sparse = TRUE)
#  #A <- as(A, "dgCMatrix")
#  # colnames(A) <- paste("q", 1:100)
#  # rownames(A) <- paste("q", 1:100)
#  dt <- data.frame(species = paste("s", mod$species, sep = ""),
#                   x = mod$x)
#  
#  Almer(formula = x ~ 1 + (1|species), data = dt, A = list(species = A))
#  A <- A[100-0:99, 100-0:99]
#  Almer(formula = x ~ 1 + (1|species), data = dt, A = list(species = A))
#  
#  mod <- Almer(formula = x ~ 1 + (1|species), data = dt, A = list(species = A))
#  
#  
#  ##### residual rate simulation #####
#  b_true <- (1:100)/20
#  b_est <- NA
#  for(i in 1:length(b_true)){
#    n <- 200
#    tree <- geiger::sim.bdtree(b = 1, d = 0, n = n, t = 4)
#    tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1]
#  
#    x <- exp(rnorm(n))
#  
#    dt <- rate_sim(tree, a = 1, b = b_true[i], x = x, model = "residual_rate")
#    mod <- rate_gls(x=dt$x, y=dt$y, species=dt$species, tree, model = "residual_rate", maxiter = 100, silent = FALSE,
#                        startv = list(a=NULL, b=NULL))
#    if(mod$convergence == "Convergence") b_est[i] <- mod$param["b", 1]
#    else b_est[i] <- NA
#  
#  }
#  
#  plot(b_true, b_est)
#  abline(0,1)
#  abline(lm(b_est~b_true[1:length(b_est)]))
#  
#  
#  mod$param
#  
#  boot.rate_gls(mod)
#  

