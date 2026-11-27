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
  expect_error(evolvabilityBeta(matrix(1, 2, 3), c(1, 0)), "not symmetric")
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
