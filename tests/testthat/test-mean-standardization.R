G <- matrix(c(4, 1, 1, 9), nrow = 2)
G2 <- matrix(c(2, 0.5, 0.5, 5), nrow = 2)
beta <- c(0.6, 0.8)

test_that("means of one are equivalent to no mean standardization", {
  expect_equal(evolvabilityMeans(G), evolvabilityMeans(G, c(1, 1)))
})

test_that("mean standardization is equivalent to standardizing G first", {
  for (means in list(c(1, 2), c(2, 1))) {
    standardized <- meanStdG(G, means)
    expect_equal(evolvabilityMeans(G, means), evolvabilityMeans(standardized))

    actual_beta <- evolvabilityBeta(G, beta, means)
    expected_beta <- evolvabilityBeta(standardized, beta)
    expect_equal(actual_beta$e, expected_beta$e)
    expect_equal(actual_beta$r, expected_beta$r)
    expect_equal(actual_beta$c, expected_beta$c)

    expect_equal(responsediffBeta(G, G2, beta, means1 = means)$d,
                 responsediffBeta(standardized, G2, beta)$d)
    expect_equal(responsediffBeta(G2, G, beta, means2 = means)$d,
                 responsediffBeta(G2, standardized, beta)$d)
  }
})

test_that("means must match the dimension of G", {
  expect_error(evolvabilityBeta(G, beta, means = c(1, 2, 3)), "length")
  expect_error(evolvabilityMeans(G, means = 2), "length")
})

test_that("meanStdGMCMC standardizes each draw with its own means", {
  G_mcmc <- rbind(c(G), c(G), c(G2))
  means_mcmc <- rbind(c(1, 1), c(1, 2), c(2, 1))
  expected_mcmc <- rbind(c(4, 1, 1, 9),
                         c(4, 0.5, 0.5, 2.25),
                         c(0.5, 0.25, 0.25, 5))
  expect_equal(meanStdGMCMC(G_mcmc, means_mcmc), expected_mcmc)
  expect_equal(meanStdGMCMC(G_mcmc[2, , drop = FALSE],
                            means_mcmc[2, , drop = FALSE]),
               expected_mcmc[2, , drop = FALSE])
})
