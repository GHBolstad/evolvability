## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE 
)
library(evolvability)

## -----------------------------------------------------------------------------
set.seed(102)
tree <- ape::rtree(n = 50)
tree <- ape::chronopl(tree, lambda = 1, age.min = 2)

## ---- eval=TRUE---------------------------------------------------------------
ape::is.ultrametric(tree)

## ---- eval=TRUE---------------------------------------------------------------
tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1]

## ---- eval = TRUE-------------------------------------------------------------
plot(tree)

## -----------------------------------------------------------------------------
sim_data <- simulate_rate(tree, startv_x=0, sigma_x=0.25, a=2, b=1, model = "predictor_BM")
head(sim_data)

## -----------------------------------------------------------------------------
gls_mod <- rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, 
                    model = "predictor_BM", silent = TRUE)
gls_mod$param

## -----------------------------------------------------------------------------
bootout <- rate_gls_boot(gls_mod, n = 10, silent = TRUE)

## -----------------------------------------------------------------------------
bootout$summary

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod, scale = "VAR")

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod) # with the default: scale == "SD"

## -----------------------------------------------------------------------------
sim_data <- simulate_rate(tree, startv_x=1, sigma_x=1, a=1, b=0.1, model = "predictor_gBM")
head(sim_data)

## -----------------------------------------------------------------------------
gls_mod <- rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "predictor_gBM", silent = TRUE)
gls_mod$param

## -----------------------------------------------------------------------------
bootout <- rate_gls_boot(gls_mod, n = 10, silent = TRUE) 

## -----------------------------------------------------------------------------
bootout$summary

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod, scale = "VAR")

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod) # with the default scale == "SD"

## -----------------------------------------------------------------------------
sim_data <- simulate_rate(tree, startv_x=0, sigma_x=1, a=1, b=1, sigma_y = 1, model = "recent_evol")
head(sim_data)

## -----------------------------------------------------------------------------
gls_mod <- rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "recent_evol", useLFO = FALSE, silent = TRUE)
gls_mod$param

## -----------------------------------------------------------------------------
bootout <- rate_gls_boot(gls_mod, n = 10, useLFO = FALSE, silent = TRUE) 

## ---- eval = TRUE-------------------------------------------------------------
bootout$summary

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod, scale = "VAR")

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod) # with the default scale == "SD"

