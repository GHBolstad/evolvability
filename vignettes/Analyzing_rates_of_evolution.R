## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE 
)
library(evolvability)

## -----------------------------------------------------------------------------
set.seed(200)
tree <- ape::rtree(n = 200)
tree <- ape::chronopl(tree, lambda = 1, age.min = 2)

## ---- eval=TRUE---------------------------------------------------------------
ape::is.ultrametric(tree)

## ---- eval=TRUE---------------------------------------------------------------
tree$edge.length <- tree$edge.length/diag(ape::vcv(tree))[1]

## ---- eval = TRUE-------------------------------------------------------------
plot(tree)

## -----------------------------------------------------------------------------
sim_data <- simulate_rate(tree, startv_x=0, sigma_x=0.25, a=1, b=1, model = "predictor_BM")

## ---- fig.width=7-------------------------------------------------------------
par(mfrow = c(1,3))
plot(sim_data) # defaults to response = "rate_y" 
plot(sim_data, response = "y")
plot(sim_data, response = "x")

## -----------------------------------------------------------------------------
d <- sim_data$tips
head(d)

## -----------------------------------------------------------------------------
gls_mod <- rate_gls(x=d$x, y=d$y, species=d$species, tree, 
                    model = "predictor_BM", silent = TRUE)
round(gls_mod$param, 3)

## -----------------------------------------------------------------------------
bootout <- rate_gls_boot(gls_mod, n = 10, silent = TRUE) 
#NB! n must be much larger to provide reliable estimates

## -----------------------------------------------------------------------------
round(bootout$summary, 3)

## ---- fig.height=5, fig.width=5-----------------------------------------------
plot(gls_mod, scale = "VAR", ylab = "y^2", xlab = "x*")

## ---- fig.height=5, fig.width=5-----------------------------------------------
plot(gls_mod, ylab = "|y|", xlab = "x*") # with the default: scale = "SD"

## ---- fig.height=5, fig.width=5-----------------------------------------------
plot(gls_mod$data$x_original, gls_mod$data$x, xlab = "x", ylab = "x*")

## -----------------------------------------------------------------------------
sim_data <- simulate_rate(tree, startv_x=1, sigma_x=1, a=1, b=0.1, model = "predictor_gBM")

## ---- fig.width=7-------------------------------------------------------------
par(mfrow = c(1,3))
plot(sim_data) # defaults to response = "rate_y" 
plot(sim_data, response = "y")
plot(sim_data, response = "x")

## -----------------------------------------------------------------------------
d <- sim_data$tips
head(d)

## -----------------------------------------------------------------------------
gls_mod <- rate_gls(x=d$x, y=d$y, species=d$species, tree, 
                    model = "predictor_gBM", silent = TRUE)
round(gls_mod$param, 3)

## -----------------------------------------------------------------------------
bootout <- rate_gls_boot(gls_mod, n = 10, silent = TRUE)
#NB! n must be much larger to provide reliable estimates
round(bootout$summary, 3)

## ---- fig.height=5, fig.width=5-----------------------------------------------
plot(gls_mod, scale = "VAR", ylab = "y^2", xlab = "x*")

## ---- fig.height=5, fig.width=5-----------------------------------------------
plot(gls_mod, ylab = "|y|", xlab = "x*") # with the default: scale = "SD"

## ---- fig.height=5, fig.width=5-----------------------------------------------
plot(gls_mod$data$x_original, gls_mod$data$x, xlab = "x", ylab = "x*")

## -----------------------------------------------------------------------------
sim_data <- simulate_rate(tree, startv_x=0, sigma_x=0.25, a=1, b=1, sigma_y = 1, model = "recent_evol")

## -----------------------------------------------------------------------------
d <- sim_data$tips
head(d)

## -----------------------------------------------------------------------------
gls_mod <- rate_gls(x=d$x, y=d$y, species=d$species, 
                    tree, model = "recent_evol", useLFO = TRUE, silent = FALSE)
round(gls_mod$param, 3)

## -----------------------------------------------------------------------------
bootout <- rate_gls_boot(gls_mod, n = 3, useLFO = FALSE, silent = TRUE) 
round(bootout$summary, 3)

## ---- fig.height=5, fig.width=5-----------------------------------------------
plot(gls_mod, scale = "VAR", ylab = "y^2", xlab = "x*")

## ---- fig.height=5, fig.width=5-----------------------------------------------
plot(gls_mod, ylab = "|y|", xlab = "x*") # with the default scale == "SD"

## ---- fig.height=5, fig.width=5-----------------------------------------------
a <- gls_mod$param["a",1]
b <- gls_mod$param["b",1]
V_micro <- a*diag(nrow(d)) + diag(b*d$x)
diag(V_micro)[diag(V_micro) < 0] <- 0  # Negative variances are replaced by zero
sigma2_y <- gls_mod$param["sigma(y)^2",1]
V_macro <- ape::vcv(tree)*sigma2_y
y_macro_pred <- macro_pred(d$y, V=V_macro+V_micro)
plot(d$y, y_macro_pred, xlab = "y", ylab = "Macro-evolutionary predictions of y")

