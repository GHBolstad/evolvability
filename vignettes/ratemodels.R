## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

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
#  rate_gls(x=mod$x, y=mod$y, species=tree$tip.label, tree, model = "predictor_geometricBM", maxiter = 100, silent = FALSE)
#  
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
#  ##### phylogenetic heritability #####
#  
#  VarCorr(mod)["Residual",]
#  
#  

## ---- eval = FALSE-------------------------------------------------------
#  ## making a dataset
#  n_species <- 10
#  tree <- geiger::sim.bdtree(b = 1, d = 0, n = n_species, t = 4)
#  tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1]
#  mod <- rate_sim(tree, startv_x=3, sigma_x=1, a= 0.001, b=1) # simulating a BM process stored in mod$x
#  mod$x <- mod$x + rnorm(n_species, 0, 0.2) # adding some residual variation
#  
#  A <- Matrix::Matrix(ape::vcv(tree), sparse = TRUE)
#  #A <- as(A, "dgCMatrix")
#  # colnames(A) <- paste("q", 1:100)
#  # rownames(A) <- paste("q", 1:100)
#  dt <- data.frame(species = paste("s", mod$species, sep = ""),
#                   x = mod$x)
#  
#  # running model
#  mod1 <- Almer(x ~ 1 + (1|species), data = dt, A = list(species = A))

## ---- eval = FALSE-------------------------------------------------------
#  # Vector of SE
#  dt$SE <- exp(rnorm(n_species))
#  dt$x <- rnorm(10, mod$x, dt$SE)
#  
#  requre(lme4) #for some reason the function does not manage to find the VarCorr function in the lme4 package
#  
#  mod_Almer <- Almer_SE(x ~ 1 + (1|species), data = dt, A = list(species = A), SE = SE)
#  
#  
#  mod_lmer <- lmer(x ~ 1 + (1|species), data = dt, control = lmerControl(check.nobs.vs.nlev  = "ignore",
#                                          check.nobs.vs.rankZ = "ignore",
#                                          check.nobs.vs.nRE   = "ignore"))
#  

