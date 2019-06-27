## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
n_species <- 100
tree <- geiger::sim.bdtree(b = 1, d = 0, n = n_species, t = 4)
tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1] 

## ------------------------------------------------------------------------
A <- Matrix::Matrix(ape::vcv(tree), sparse = TRUE)

## ------------------------------------------------------------------------
colnames(A) <- rownames(A) <- paste("species", 1:n_species, sep = "_")

## ------------------------------------------------------------------------
y <- 5 + t(chol(A))%*%rnorm(n_species, 0, 2) + # BM process with mean 5 and s.d. 2
     rnorm(n_species, 0, 1)             # residual variation with s.d. 1

## ------------------------------------------------------------------------
dt <- data.frame(species = colnames(A),
                 y = as.vector(y))

## ------------------------------------------------------------------------
mod <- Almer(y ~ 1 + (1|species), data = dt, A = list(species = A))
mod

## ---- eval = FALSE-------------------------------------------------------
#  dt$SE <- runif(nrow(dt), min = 0.5, max = 1) # adding some arbitrary SE values
#  
#  mod_SE <- Almer_SE(y ~ 1 + (1|species), data = dt, SE = dt$SE, A = list(species = A))
#  mod_SE

## ------------------------------------------------------------------------
sim_y <- simulate_Almer(mod, nsim = 3)
head(sim_y)

## ------------------------------------------------------------------------
phylH(mod, numerator = "species", nsim = 10)

