asymmetric <- matrix(c(1, 5, 0, 1, 2, 2, 0, 2, 3), nrow = 3)
beta <- c(1, 0, 0)

test_that("asymmetric G matrices are rejected", {
  expect_error(evolvabilityBeta(asymmetric, beta), "G is not symmetric")
  expect_error(evolvabilityMeans(asymmetric), "G is not symmetric")
  expect_error(conditionalG(asymmetric, 1), "G is not symmetric")
  expect_error(responsediffBeta(asymmetric, diag(3), beta),
               "G1 is not symmetric")
  expect_error(responsediffBeta(diag(3), asymmetric, beta),
               "G2 is not symmetric")
  expect_error(evolvabilityBeta(matrix(1, 2, 3), c(1, 0)), "square matrix")
  expect_error(evolvabilityMeans(c(1, 2, 2, 1)), "square matrix")
})

test_that("rounding errors in G are tolerated", {
  G <- crossprod(matrix(c(1, 2, 0, 1, 3, 1, 0, 1, 2), 3))
  G[1, 2] <- G[1, 2] + 1e-12
  expect_no_error(evolvabilityBeta(G, beta))
})

test_that("conditionalG conditions on the selected traits", {
  G <- matrix(c(4, 1, 2, 1, 3, 0.5, 2, 0.5, 5), nrow = 3)
  expected <- G[2:3, 2:3] - G[2:3, 1, drop = FALSE] %*% t(G[2:3, 1, drop = FALSE]) / G[1, 1]
  expect_equal(conditionalG(G, 1), expected)
  # a single remaining trait is returned as a 1 x 1 matrix
  single <- conditionalG(G, c(1, 3))
  expect_equal(dim(single), c(1, 1))
  expect_equal(single[1, 1], 1 / solve(G)[2, 2])
  expect_error(conditionalG(G), "Specify")
  expect_error(conditionalG(G, 4), "column numbers")
  expect_error(conditionalG(G, 1:3), "all the columns")
})

test_that("non-finite and malformed G matrices give clear errors", {
  with_na <- diag(2)
  with_na[1, 1] <- NA
  expect_error(evolvabilityBeta(with_na, c(1, 0)), "missing or non-finite")
  expect_error(evolvabilityMeans(with_na), "missing or non-finite")
  expect_error(conditionalG(with_na, 1), "missing or non-finite")
  expect_error(evolvabilityBeta(diag(3), c(1, 0)), "as many rows")
})

test_that("MCMC functions reject asymmetric or malformed G draws", {
  G <- diag(c(2, 3))
  draws <- matrix(rep(c(G), 5), nrow = 5, byrow = TRUE)
  asymmetric <- draws
  asymmetric[4, ] <- c(2, 1, 0, 3)
  beta <- c(1, 0)
  beta_draws <- matrix(rep(beta, 5), nrow = 5, byrow = TRUE)

  expect_error(evolvabilityBetaMCMC(asymmetric, beta), "not symmetric.*row 4")
  expect_error(evolvabilityBetaMCMC2(asymmetric, beta_draws),
               "not symmetric.*row 4")
  expect_error(evolvabilityMeansMCMC(asymmetric), "not symmetric.*row 4")

  expect_error(evolvabilityBetaMCMC(draws[, 1:3], beta), "square")
  expect_error(evolvabilityBetaMCMC(draws, c(1, 0, 0)), "as many rows")
  expect_error(evolvabilityBetaMCMC2(draws, beta_draws[1:3, ]),
               "one row per row")
  draws[2, 1] <- NA
  expect_error(evolvabilityBetaMCMC(draws, beta), "missing or non-finite")
})

test_that("conditionalG returns NA with a warning for a singular conditioning block", {
  expect_warning(result <- conditionalG(diag(c(2, 0, 1)), 2),
                 "could not be inverted by solve()", fixed = TRUE)
  expect_equal(dim(result), c(2, 2))
  expect_true(all(is.na(result)))
  # conditioning on a non-singular block of a singular G still works
  expect_no_warning(result <- conditionalG(diag(c(2, 0, 1)), 1))
  expect_equal(result, diag(c(0, 1)))
})

test_that("meanStdG and meanStdGMCMC reject asymmetric G", {
  expect_error(meanStdG(asymmetric, c(1, 2, 3)), "G is not symmetric")
  expect_error(meanStdG(diag(3), c(1, 2)), "length of the means")
  draws <- rbind(c(diag(2)), c(2, 1, 0, 2))
  means <- rbind(c(1, 1), c(1, 1))
  expect_error(meanStdGMCMC(draws, means), "not symmetric.*row 2")
  expect_error(meanStdGMCMC(rbind(c(diag(2)), c(diag(2))), means[, 1, drop = FALSE]),
               "means_mcmc")
})
