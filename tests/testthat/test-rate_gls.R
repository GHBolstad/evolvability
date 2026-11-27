skip_on_cran()

test_that("rate_gls fits the predictor_BM model and reports progress", {
  set.seed(4)
  tree <- ape::rcoal(40)
  tree$edge.length <- tree$edge.length / max(ape::node.depth.edgelength(tree))
  sim <- simulate_rate(tree, startv_x = 0, sigma_x = 0.25, a = 1, b = 1,
                       model = "predictor_BM")$tips

  fit <- rate_gls(sim$x, sim$y, sim$species, tree, silent = TRUE)
  expect_s3_class(fit, "rate_gls")
  expect_equal(fit$convergence, "Convergence")
  expect_equal(rownames(fit$param), c("a", "b", "sigma(x)^2"))

  expect_message(rate_gls(sim$x, sim$y, sim$species, tree, silent = FALSE),
                 "converg score")
  expect_silent(rate_gls(sim$x, sim$y, sim$species, tree, silent = TRUE))
})
