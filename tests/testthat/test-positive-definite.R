# A symmetric matrix that is invertible but has a negative eigenvalue
not_pd <- matrix(c(1, 2, 2, 1), nrow = 2)
pd <- diag(c(2, 3))
beta <- c(0.6, 0.8)

expect_not_pd_warning <- function(expr) {
  value <- NULL
  expect_warning(value <- expr, "not positive definite")
  value
}

test_that("evolvabilityBeta warns about non-positive definite G but calculates", {
  result <- expect_not_pd_warning(evolvabilityBeta(not_pd, beta))
  expect_true(all(is.finite(c(result$e, result$r, result$c, result$a, result$i))))
  expect_no_warning(evolvabilityBeta(pd, beta))
})

test_that("a singular G that is positive semi-definite only gives the inversion warning", {
  warnings <- character()
  withCallingHandlers(
    evolvabilityBeta(diag(c(2, 0)), beta),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings, 1)
  expect_match(warnings, "could not be inverted by solve()", fixed = TRUE)
})

test_that("evolvabilityMeans warns about non-positive definite G but calculates", {
  result <- expect_not_pd_warning(evolvabilityMeans(not_pd))
  expect_true(all(is.finite(result)))
  expect_equal(unname(result[c("e_min", "e_max")]), c(-1, 3))
})

test_that("conditionalG and responsediffBeta warn about non-positive definite G", {
  expect_not_pd_warning(conditionalG(diag(3) + matrix(c(0, 2, 0, 2, 0, 0, 0, 0, 0), 3), 3))
  result <- NULL
  expect_warning(result <- responsediffBeta(not_pd, pd, beta),
                 "G1 is not positive definite")
  expect_equal(unname(result$d),
               sqrt(sum(((not_pd - pd) %*% beta)^2)))
  expect_warning(responsediffBeta(pd, not_pd, beta),
                 "G2 is not positive definite")
})

test_that("the MCMC functions warn about non-positive definite draws", {
  draws <- matrix(rep(c(pd), 10), nrow = 10, byrow = TRUE)
  draws[4, ] <- c(not_pd)
  beta_draws <- matrix(rep(beta, 10), nrow = 10, byrow = TRUE)

  result <- expect_not_pd_warning(evolvabilityBetaMCMC(draws, beta, post.dist = TRUE))
  expect_true(all(is.finite(result$post.dist$cB)))
  result <- expect_not_pd_warning(evolvabilityBetaMCMC2(draws, beta_draws,
                                                        post.dist = TRUE))
  expect_true(all(is.finite(result$post.dist)))
  result <- expect_not_pd_warning(evolvabilityMeansMCMC(draws))
  expect_true(all(is.finite(result$post.dist)))
})

test_that("conditional parameters are NA and others still calculated in every function", {
  singular_g <- diag(c(2, 0))
  draws <- matrix(rep(c(singular_g), 5), nrow = 5, byrow = TRUE)
  beta_draws <- matrix(rep(beta, 5), nrow = 5, byrow = TRUE)
  inverted <- "could not be inverted by solve()"

  result <- suppressWarnings(evolvabilityBeta(singular_g, beta))
  expect_true(is.finite(result$e) && is.finite(result$r))
  expect_true(is.na(result$c) && is.na(result$a) && is.na(result$i))

  result <- suppressWarnings(evolvabilityBetaMCMC(draws, beta))
  expect_true(all(is.finite(c(result$eB, result$rB))))
  expect_true(all(is.na(c(result$cB, result$aB, result$iB))))

  result <- suppressWarnings(evolvabilityBetaMCMC2(draws, beta_draws))
  expect_true(all(is.finite(result$summary[c("eB", "rB"), ])))
  expect_true(all(is.na(result$summary[c("cB", "aB", "iB"), ])))

  expect_warning(evolvabilityBetaMCMC2(draws, beta_draws), inverted, fixed = TRUE)
  expect_warning(evolvabilityMeansMCMC(draws), inverted, fixed = TRUE)
})
