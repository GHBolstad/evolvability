# Expects exactly the inversion warning and returns the value of expr
expect_inversion_warning <- function(expr) {
  value <- NULL
  expect_warning(value <- expr, "could not be inverted by solve()", fixed = TRUE)
  value
}

G <- diag(c(2, 3))
singular <- diag(c(2, 0))
near_singular <- diag(c(1, 1e-20))
Beta <- diag(2)

G_mcmc <- matrix(rep(c(G), 10), nrow = 10, byrow = TRUE)
G_mcmc[3, ] <- c(singular)
draws <- matrix(rep(c(1, 0), 10), nrow = 10, byrow = TRUE)

test_that("evolvabilityBeta gives NA for conditional parameters of singular G", {
  ordinary <- evolvabilityBeta(G, Beta)
  expect_equal(unname(ordinary$c), c(2, 3))

  for (bad in list(singular, near_singular)) {
    result <- expect_inversion_warning(evolvabilityBeta(bad, Beta))
    expect_equal(unname(result$e), diag(bad))
    expect_equal(unname(result$r), diag(bad))
    expect_true(all(is.na(result$c)))
    expect_true(all(is.na(result$a)))
    expect_true(all(is.na(result$i)))
    expect_true(all(is.na(summary(result)$Averages[3:5])))
  }
})

test_that("evolvabilityBetaMCMC handles singular draws", {
  result <- expect_inversion_warning(
    evolvabilityBetaMCMC(G_mcmc, Beta, post.dist = TRUE)
  )
  expect_true(all(is.finite(result$post.dist$eB)))
  expect_true(all(is.finite(result$post.dist$rB)))
  expect_true(all(is.na(result$post.dist$cB[3, ])))
  expect_true(all(is.finite(result$post.dist$cB[-3, ])))
  expect_true(all(is.na(result$cB)))
  expect_true(all(is.na(result$aB)))
  expect_true(all(is.na(result$iB)))
  expect_true(all(is.finite(result$eB)))
  expect_true(all(is.finite(result$rB)))
  expect_true(all(is.na(result$summary[, 3:5])))
  expect_true(all(is.finite(result$summary[, 1:2])))

  summarized <- summary(result)
  expect_true(all(is.na(summarized$Minimum[, 3:5])))
  expect_true(all(is.na(summarized$Maximum[, 3:5])))
})

test_that("evolvabilityBetaMCMC2 handles singular draws", {
  result <- expect_inversion_warning(
    evolvabilityBetaMCMC2(G_mcmc, draws, post.dist = TRUE)
  )
  expect_true(all(is.finite(result$post.dist[, 1:2])))
  expect_true(all(is.na(result$post.dist[3, 3:5])))
  expect_true(all(is.finite(result$post.dist[-3, 3:5])))
  expect_true(all(is.na(result$summary[3:5, ])))
  expect_true(all(is.finite(result$summary[1:2, ])))
})

test_that("all draws singular", {
  G_all <- matrix(rep(c(singular), 10), nrow = 10, byrow = TRUE)
  result <- expect_inversion_warning(
    evolvabilityBetaMCMC(G_all, c(1, 0), post.dist = TRUE)
  )
  expect_true(all(is.na(result$post.dist$cB)))
  expect_true(all(is.finite(result$post.dist$eB)))

  result <- expect_inversion_warning(
    evolvabilityBetaMCMC(G_all, cbind(c(1, 0), c(1, 0)))
  )
  expect_true(all(is.na(summary(result)$Minimum[, 3:5])))

  result <- expect_inversion_warning(
    evolvabilityBetaMCMC2(G_all, draws, post.dist = TRUE)
  )
  expect_true(all(is.na(result$post.dist[, 3:5])))
  expect_true(all(is.finite(result$post.dist[, 1:2])))
})

test_that("solve_g only treats singularity, not bad input, as a failed inversion", {
  solve_g <- evolvability:::solve_g
  expect_null(solve_g(singular))
  expect_equal(solve_g(G), diag(1 / c(2, 3)))
  expect_error(solve_g(matrix(1:6, 2, 3)))
})

test_that("evolvabilityMeans warns for singular G and returns what it can", {
  for (bad in list(singular, near_singular)) {
    result <- expect_inversion_warning(evolvabilityMeans(bad))
    expect_true(all(is.na(result[c("c_mean", "a_mean", "i_mean")])))
    expect_equal(unname(result[c("e_mean", "e_min", "e_max")]),
                 c(mean(diag(bad)), min(diag(bad)), max(diag(bad))))
    expect_true(is.finite(result[["r_mean"]]))
  }
  expect_no_warning(evolvabilityMeans(G))
})

test_that("evolvabilityMeansMCMC handles singular draws", {
  result <- expect_inversion_warning(evolvabilityMeansMCMC(G_mcmc))
  expect_equal(nrow(result$post.medians), 7)
  expect_true(all(is.na(result$post.medians[c("c_mean", "a_mean", "i_mean"), ])))
  expect_true(all(is.finite(result$post.medians[c("e_mean", "e_min", "e_max",
                                                  "r_mean"), ])))
  expect_true(all(is.na(result$post.dist[3, c("c_mean", "a_mean", "i_mean")])))
  expect_true(all(is.finite(result$post.dist[-3, ])))

  all_singular <- matrix(rep(c(singular), 10), nrow = 10, byrow = TRUE)
  result <- expect_inversion_warning(evolvabilityMeansMCMC(all_singular))
  expect_true(all(is.finite(result$post.medians[c("e_mean", "r_mean"), ])))
})

test_that("evolvabilityMeansMCMC is unchanged when all draws are invertible", {
  draws <- rbind(c(diag(c(2, 3))), c(diag(c(1, 4))), c(diag(c(3, 1))))
  res <- evolvabilityMeansMCMC(draws)
  expect_equal(unname(res$post.medians[, "median"]),
               unname(apply(res$post.dist, 2, median)))
  expect_equal(unname(res$post.medians[, c("lower", "upper")]),
               unname(coda::HPDinterval(res$post.dist)),
               ignore_attr = TRUE)
})
