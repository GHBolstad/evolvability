## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE #NB! må kjøre modellane på nytt pga endring til A og B i staden for intercept og slope
)
library(evolvability)

## -----------------------------------------------------------------------------
#  tree <- ape::rtree(n = 10)
#  tree <- ape::chronopl(tree, lambda = 1, age.min = 2)

## ---- echo=FALSE--------------------------------------------------------------
#  save(tree, file = "../data/tree.RData") # the model output is stored as it takes too long time to run within the vignette

## ---- echo=FALSE, eval = TRUE-------------------------------------------------
load("../data/tree.RData")

## ---- eval = TRUE-------------------------------------------------------------
ape::is.ultrametric(tree)

## -----------------------------------------------------------------------------
#  tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1]

## -----------------------------------------------------------------------------
#  sim_data <- simulate_rate(tree, startv_x=0, sigma_x=0.25, a=1, b=1, model = "predictor_BM")

## -----------------------------------------------------------------------------
#  gls_mod <- rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "predictor_BM")

## ---- echo=FALSE--------------------------------------------------------------
#  save(gls_mod, file = "../data/mod_BM.RData") # the model output is stored as it takes too long time to run within the vignette

## ---- echo=FALSE, eval = TRUE-------------------------------------------------
load("../data/mod_BM.RData")

## ---- eval = TRUE-------------------------------------------------------------
gls_mod$param

## -----------------------------------------------------------------------------
#  bootout <- boot.rate_gls(gls_mod, n = 10)

## ---- echo=FALSE--------------------------------------------------------------
#  save(bootout, file = "../data/mod_BM_boot.RData") # the bootstrap output is stored as it takes too long time to run within the vignette

## ---- echo=FALSE, eval = TRUE-------------------------------------------------
load("../data/mod_BM_boot.RData")

## ---- eval = TRUE-------------------------------------------------------------
bootout$summary

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod, scale = "VAR")

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod) # with the default: scale == "SD"

## -----------------------------------------------------------------------------
#  sim_data <- simulate_rate(tree, startv_x=1, sigma_x=2, a=2, b=0.1, model = "predictor_gBM")

## -----------------------------------------------------------------------------
#  gls_mod <- rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "predictor_gBM")

## ---- echo=FALSE--------------------------------------------------------------
#  save(gls_mod, file = "../data/mod_gBM.RData") # the model output is stored as it takes too long time to run within the vignette

## ---- echo=FALSE, eval = TRUE-------------------------------------------------
load("../data/mod_gBM.RData")

## ---- eval = TRUE-------------------------------------------------------------
gls_mod$param

## -----------------------------------------------------------------------------
#  bootout <- boot.rate_gls(gls_mod, n = 10)

## ---- echo=FALSE--------------------------------------------------------------
#  save(bootout, file = "../data/mod_gBM_boot.RData") # the bootstrap output is stored as it takes too long time to run within the vignette

## ---- echo=FALSE, eval = TRUE-------------------------------------------------
load("../data/mod_gBM_boot.RData")

## ---- eval = TRUE-------------------------------------------------------------
bootout$summary

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod, scale = "VAR")

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod) # with the default scale == "SD"

## -----------------------------------------------------------------------------
#  sim_data <- simulate_rate(tree, startv_x=0, sigma_x=1, a=1, b=1, sigma_y = 1, model = "recent_evol")

## -----------------------------------------------------------------------------
#  gls_mod <- rate_gls(x=sim_data$x, y=sim_data$y, species=sim_data$species, tree, model = "recent_evol", useLFO = FALSE)

## ---- echo=FALSE--------------------------------------------------------------
#  save(gls_mod, file = "../data/mod_recentevol.RData") # the model output is stored as it takes too long time to run within the vignette

## ---- echo=FALSE, eval = TRUE-------------------------------------------------
load("../data/mod_recentevol.RData")

## ---- eval = TRUE-------------------------------------------------------------
gls_mod$param

## -----------------------------------------------------------------------------
#  bootout <- boot.rate_gls(gls_mod, n = 10)

## ---- echo=FALSE--------------------------------------------------------------
#  save(bootout, file = "../data/mod_recentevol_boot.RData") # the bootstrap output is stored as it takes too long time to run within the vignette

## ---- echo=FALSE, eval = TRUE-------------------------------------------------
load("../data/mod_recentevol_boot.RData")

## ---- eval = TRUE-------------------------------------------------------------
bootout$summary

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod, scale = "VAR")

## ---- fig.height=5, fig.width=5, eval = TRUE----------------------------------
plot(gls_mod) # with the default scale == "SD"

