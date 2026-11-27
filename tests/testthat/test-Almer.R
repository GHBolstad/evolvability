dat <- lme4::sleepstudy
form <- Reaction ~ Days + (1 | Subject)

test_that("weights and offset are forwarded to lme4", {
  # The formula environment must contain the weights and the offset
  form <- Reaction ~ Days + (1 | Subject)
  w <- seq(1, 2, length.out = nrow(dat))
  off <- dat$Days / 2
  unweighted <- Almer(form, data = dat)
  weighted <- Almer(form, data = dat, weights = w, offset = off)
  reference <- lme4::lmer(form, data = dat, weights = w, offset = off)

  expect_null(stats::model.weights(stats::model.frame(unweighted)))
  expect_identical(stats::model.weights(stats::model.frame(weighted)), w)
  expect_identical(stats::model.offset(stats::model.frame(weighted)), off)
  expect_equal(lme4::fixef(weighted), lme4::fixef(reference), tolerance = 1e-5)
})

test_that("variables are found in the calling environment", {
  unweighted <- Almer(form, data = dat)
  local_fit <- local({
    response <- dat$Reaction
    predictor <- dat$Days
    group <- dat$Subject
    Almer(response ~ predictor + (1 | group))
  })
  expect_equal(lme4::fixef(local_fit), lme4::fixef(unweighted),
               ignore_attr = TRUE)
})

test_that("an A matrix that matches no random effect gives a clear error", {
  expect_error(
    Almer(form, data = dat, A = list(Day = diag(nlevels(dat$Subject)))),
    "does not match any random effect"
  )
})

test_that("Almer_SE warns only when it fails to converge", {
  SE <- rep(5, nrow(dat))
  expect_warning(Almer_SE(form, SE = SE, data = dat, maxiter = 1),
                 "maximum number of iterations")
  expect_no_warning(Almer_SE(form, SE = SE, data = dat, maxiter = 100))
})

test_that("Almer_boot and phylH return bootstrap summaries", {
  mod <- Almer(form, data = dat)
  set.seed(1)
  boot <- Almer_boot(mod, nsim = 3)
  expect_equal(dim(boot$fixef_distribution), c(3, 2))
  expect_equal(dim(boot$vcov_distribution), c(3, 2))

  set.seed(1)
  h <- phylH(mod, "Subject", nsim = 3)
  expect_named(h$bootstrap, paste0("Sim_", 1:3))
  expect_length(h$phylH, 3)
})
