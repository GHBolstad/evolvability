## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
n_species <- 100
tree <- ape::rtree(n = n_species)
tree <- ape::chronopl(tree, lambda = 1)

## -----------------------------------------------------------------------------
A <- Matrix::Matrix(ape::vcv(tree), sparse = TRUE)

## -----------------------------------------------------------------------------
colnames(A) <- rownames(A) <- paste("species", 1:n_species, sep = "_")

## -----------------------------------------------------------------------------
y <- 5 + t(chol(A))%*%rnorm(n_species, 0, 2) + # BM process with mean 5 and s.d. 2
     rnorm(n_species, 0, 1)                    # residual variation with s.d. 1

## -----------------------------------------------------------------------------
dt <- data.frame(species = colnames(A),
                 y = as.vector(y))

## -----------------------------------------------------------------------------
mod <- Almer(y ~ 1 + (1|species), data = dt, A = list(species = A))
summary(mod)

## -----------------------------------------------------------------------------
dt$SE <- runif(nrow(dt), min = 0.1, max = 0.2) 

## -----------------------------------------------------------------------------
mod_SE <- Almer_SE(y ~ 1 + (1|species), data = dt, SE = dt$SE, A = list(species = A))
summary(mod_SE)

## -----------------------------------------------------------------------------
sim_y <- Almer_sim(mod_SE, nsim = 3)
sim_y[1:3,]

## ---- eval = FALSE------------------------------------------------------------
#  dt$group <- c(rep("A", 50), rep("B", 50))
#  mod <- Almer(y ~ group + (1|species), data = dt, A = list(species = A))
#  #mod <- Almer(y ~ 1 + (1|species), data = dt, A = list(species = A))
#  nsim = 3
#  y <- Almer_sim(mod, nsim)
#  dt <- mod@frame
#  boot <- apply(y, 2, function(x){
#    dt$hopefullynotavariablealready <- x
#    update(mod, hopefullynotavariablealready ~ ., data = dt)
#  })
#  fixef_boot <- sapply(boot, lme4::fixef)
#  if(is.vector(fixef_boot)){
#    fixef_boot <- cbind(fixef_boot)
#  }else fixef_boot <- t(fixef_boot)
#  colnames(fixef_boot) <- names(lme4::fixef(boot[[1]]))
#  
#  VarCorr_list <- lapply(boot, function(x) as.data.frame(lme4::VarCorr(x)))
#  vcov_boot <- t(sapply(VarCorr_list, function(x) x$vcov))
#  colnames(vcov_boot) <- VarCorr_list[[1]]$grp
#  
#  list(
#    fixef = cbind(mean = apply(fixef_boot, 2, mean),
#        std.err. = apply(fixef_boot, 2, sd),
#        "2.5% quantile" = apply(fixef_boot, 2, quantile, probs = 0.025),
#        "97.5% quantile" = apply(fixef_boot, 2, quantile, probs = 0.975)
#        ),
#    vcov = cbind(mean = apply(vcov_boot, 2, mean),
#        std.err. = apply(vcov_boot, 2, sd),
#        "2.5% quantile" = apply(vcov_boot, 2, quantile, probs = 0.025),
#        "97.5% quantile" = apply(vcov_boot, 2, quantile, probs = 0.975)
#        ),
#    models <- boot)
#  
#  

## -----------------------------------------------------------------------------
phylH(mod, numerator = "species", nsim = 3)

