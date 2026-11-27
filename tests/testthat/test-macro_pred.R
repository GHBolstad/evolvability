# Direct implementation of the macroevolutionary predictions
manual_macro_pred <- function(y, V, lfo) {
  gls_mean <- function(y, V) {
    Vi <- solve(V)
    sum(Vi %*% y) / sum(Vi)
  }
  y_mean <- if (lfo) {
    vapply(seq_along(y), function(i) gls_mean(y[-i], V[-i, -i]), numeric(1))
  } else {
    rep(gls_mean(y, V), length(y))
  }
  Vinv <- solve(V)
  c(y_mean + (V - diag(1 / diag(Vinv))) %*% Vinv %*% (y - y_mean))
}

set.seed(10)
V <- ape::vcv(ape::rtree(6))
y <- rnorm(6)

test_that("macro_pred matches a direct implementation", {
  expect_equal(c(as.matrix(macro_pred(y, V))), manual_macro_pred(y, V, TRUE))
  expect_equal(c(as.matrix(macro_pred(y, V, useLFO = FALSE))),
               manual_macro_pred(y, V, FALSE))
})

test_that("macro_pred handles tied trait values", {
  y_tied <- y
  y_tied[2] <- y_tied[1]
  expect_equal(c(as.matrix(macro_pred(y_tied, V))),
               manual_macro_pred(y_tied, V, TRUE))
})
