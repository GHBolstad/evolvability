library(evolvability)

check_forwarded_arguments <- function() {
  dat <- lme4::sleepstudy
  form <- Reaction ~ Days + (1 | Subject)
  w <- seq(1, 2, length.out = nrow(dat))
  off <- dat$Days / 2
  unweighted <- Almer(form, data = dat)
  weighted <- Almer(form, data = dat, weights = w, offset = off)
  reference <- lme4::lmer(form, data = dat, weights = w, offset = off)
  stopifnot(
    is.null(stats::model.weights(stats::model.frame(unweighted))),
    identical(stats::model.weights(stats::model.frame(weighted)), w),
    identical(stats::model.offset(stats::model.frame(weighted)), off),
    isTRUE(all.equal(lme4::fixef(weighted), lme4::fixef(reference),
                     tolerance = 1e-5))
  )

  local_fit <- local({
    response <- dat$Reaction
    predictor <- dat$Days
    group <- dat$Subject
    Almer(response ~ predictor + (1 | group))
  })
  stopifnot(isTRUE(all.equal(lme4::fixef(local_fit),
                             lme4::fixef(unweighted),
                             check.attributes = FALSE)))
}

check_forwarded_arguments()
