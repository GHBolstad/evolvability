library(evolvability)

check_warnings <- function(expr, expected) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  stopifnot(length(warnings) == expected)
  if (expected > 0) {
    stopifnot(all(warnings == "Selection gradients (in Beta) are not of unit length"))
  }
  value
}

G <- diag(c(2, 3))
G_mcmc <- matrix(rep(c(G), 10), nrow = 10, byrow = TRUE)
unit <- c(0.6, 0.8)
near_unit <- unit * (1 + sqrt(.Machine$double.eps) / 2)
non_unit <- unit * 2
gradients <- cbind(unit, non_unit)

check_warnings(evolvabilityBeta(G, unit), 0)
check_warnings(evolvabilityBeta(G, near_unit), 0)
result <- check_warnings(evolvabilityBeta(G, gradients), 1)
stopifnot(isTRUE(all.equal(unname(result$e), c(2.64, 10.56))))
check_warnings(evolvabilityBeta(G, c(0, 0)), 1)
check_warnings(responsediffBeta(G, diag(2), unit), 0)
result <- check_warnings(responsediffBeta(G, diag(2), non_unit), 1)
stopifnot(isTRUE(all.equal(unname(result$d),
                           sqrt(sum(((G - diag(2)) %*% non_unit)^2)))))

check_warnings(evolvabilityBetaMCMC(G_mcmc, unit), 0)
result <- check_warnings(evolvabilityBetaMCMC(G_mcmc, non_unit), 1)
stopifnot(isTRUE(all.equal(unname(result$eB[1, 1]), 10.56)))
check_warnings(evolvabilityBetaMCMC(G_mcmc, gradients), 1)

draws <- matrix(rep(unit, 10), nrow = 10, byrow = TRUE)
check_warnings(evolvabilityBetaMCMC2(G_mcmc, draws), 0)
draws[5, ] <- non_unit
result <- check_warnings(evolvabilityBetaMCMC2(G_mcmc, draws, post.dist = TRUE), 1)
stopifnot(isTRUE(all.equal(unname(result$post.dist[5, "eB"]), 10.56)))

check_length <- getFromNamespace("check_beta_length", "evolvability")
check_warnings(check_length(c(NA_real_, 1)), 1)
check_warnings(check_length(c(Inf, 1)), 1)
