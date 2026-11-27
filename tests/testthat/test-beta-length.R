msg <- "Selection gradients (in Beta) are not of unit length"

# Expects the unit length warning and returns the value of expr
warns_non_unit <- function(expr) {
  value <- NULL
  expect_warning(value <- expr, msg, fixed = TRUE)
  value
}

G <- diag(c(2, 3))
G_mcmc <- matrix(rep(c(G), 10), nrow = 10, byrow = TRUE)
unit <- c(0.6, 0.8)
near_unit <- unit * (1 + sqrt(.Machine$double.eps) / 2)
non_unit <- unit * 2
gradients <- cbind(unit, non_unit)

test_that("evolvabilityBeta warns about non-unit selection gradients", {
  expect_no_warning(evolvabilityBeta(G, unit))
  expect_no_warning(evolvabilityBeta(G, near_unit))
  result <- warns_non_unit(evolvabilityBeta(G, gradients))
  expect_equal(unname(result$e), c(2.64, 10.56))
  expect_warning(evolvabilityBeta(G, c(0, 0)), msg, fixed = TRUE)
})

test_that("responsediffBeta warns about non-unit selection gradients", {
  expect_no_warning(responsediffBeta(G, diag(2), unit))
  result <- warns_non_unit(responsediffBeta(G, diag(2), non_unit))
  expect_equal(unname(result$d),
               sqrt(sum(((G - diag(2)) %*% non_unit)^2)))
})

test_that("evolvabilityBetaMCMC warns about non-unit selection gradients", {
  expect_no_warning(evolvabilityBetaMCMC(G_mcmc, unit))
  result <- warns_non_unit(evolvabilityBetaMCMC(G_mcmc, non_unit))
  expect_equal(unname(result$eB[1, 1]), 10.56)
  expect_warning(evolvabilityBetaMCMC(G_mcmc, gradients), msg, fixed = TRUE)
})

test_that("evolvabilityBetaMCMC2 warns about non-unit selection gradients", {
  draws <- matrix(rep(unit, 10), nrow = 10, byrow = TRUE)
  expect_no_warning(evolvabilityBetaMCMC2(G_mcmc, draws))
  draws[5, ] <- non_unit
  result <- warns_non_unit(
    evolvabilityBetaMCMC2(G_mcmc, draws, post.dist = TRUE)
  )
  expect_equal(unname(result$post.dist[5, "eB"]), 10.56)
})

test_that("check_beta_length treats non-finite values as non-unit", {
  expect_warning(evolvability:::check_beta_length(c(NA_real_, 1)), msg,
                 fixed = TRUE)
  expect_warning(evolvability:::check_beta_length(c(Inf, 1)), msg,
                 fixed = TRUE)
})
