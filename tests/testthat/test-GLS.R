set.seed(3)
n <- 12
X <- cbind(1, rnorm(n))
y <- rnorm(n)
R <- crossprod(matrix(rnorm(n * n), n)) / n + diag(n)

test_that("GLS reproduces the generalized least squares estimator", {
  Rinv <- solve(R)
  beta <- solve(t(X) %*% Rinv %*% X, t(X) %*% Rinv %*% y)
  fit <- GLS(y, X, R)
  expect_equal(unname(fit$coef[, "Estimate"]), c(beta))
  sigma2 <- c(t(y - X %*% beta) %*% Rinv %*% (y - X %*% beta)) / (n - 2)
  expect_equal(fit$sigma2, sigma2)
  expect_equal(unname(fit$coef[, "SE"]),
               sqrt(diag(sigma2 * solve(t(X) %*% Rinv %*% X))))
  expect_equal(GLS(y, X, R, coef_only = TRUE)$coef, fit$coef[, "Estimate"])
})

test_that("GLS with an identity R equals lm", {
  fit <- GLS(y, X, diag(n))
  reference <- summary(lm(y ~ X[, 2]))
  expect_equal(unname(fit$coef), unname(reference$coefficients[, 1:2]))
  expect_equal(fit$R2, reference$r.squared)
})

test_that("GLS accepts the Cholesky factor and requires R or L", {
  expect_equal(GLS(y, X, L = t(chol(R)))$coef, GLS(y, X, R)$coef)
  expect_error(GLS(y, X), "Provide either R or L")
})
