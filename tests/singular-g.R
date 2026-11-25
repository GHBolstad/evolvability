library(evolvability)

capture_warning <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  stopifnot(length(warnings) == 1,
            grepl("could not be inverted by solve()", warnings, fixed = TRUE))
  value
}

G <- diag(c(2, 3))
singular <- diag(c(2, 0))
near_singular <- diag(c(1, 1e-20))
Beta <- diag(2)

ordinary <- evolvabilityBeta(G, Beta)
stopifnot(isTRUE(all.equal(unname(ordinary$c), c(2, 3))))
for (bad in list(singular, near_singular)) {
  result <- capture_warning(evolvabilityBeta(bad, Beta))
  stopifnot(isTRUE(all.equal(unname(result$e), diag(bad))),
            isTRUE(all.equal(unname(result$r), diag(bad))),
            all(is.na(result$c)), all(is.na(result$a)), all(is.na(result$i)),
            all(is.na(summary(result)$Averages[3:5])))
}

G_mcmc <- matrix(rep(c(G), 10), nrow = 10, byrow = TRUE)
G_mcmc[3, ] <- c(singular)
draws <- matrix(rep(c(1, 0), 10), nrow = 10, byrow = TRUE)

result <- capture_warning(evolvabilityBetaMCMC(G_mcmc, Beta, post.dist = TRUE))
stopifnot(all(is.finite(result$post.dist$eB)),
          all(is.finite(result$post.dist$rB)),
          all(is.na(result$post.dist$cB[3, ])),
          all(is.finite(result$post.dist$cB[-3, ])),
          all(is.na(result$cB)), all(is.na(result$aB)), all(is.na(result$iB)),
          all(is.finite(result$eB)), all(is.finite(result$rB)),
          all(is.na(result$summary[, 3:5])),
          all(is.finite(result$summary[, 1:2])))
summarized <- summary(result)
stopifnot(all(is.na(summarized$Minimum[, 3:5])),
          all(is.na(summarized$Maximum[, 3:5])))

result <- capture_warning(evolvabilityBetaMCMC2(G_mcmc, draws, post.dist = TRUE))
stopifnot(all(is.finite(result$post.dist[, 1:2])),
          all(is.na(result$post.dist[3, 3:5])),
          all(is.finite(result$post.dist[-3, 3:5])),
          all(is.na(result$summary[3:5, ])),
          all(is.finite(result$summary[1:2, ])))

G_mcmc[,] <- matrix(rep(c(singular), 10), nrow = 10, byrow = TRUE)
result <- capture_warning(evolvabilityBetaMCMC(G_mcmc, c(1, 0), post.dist = TRUE))
stopifnot(all(is.na(result$post.dist$cB)), all(is.finite(result$post.dist$eB)))
result <- capture_warning(evolvabilityBetaMCMC(G_mcmc, cbind(c(1, 0), c(1, 0))))
stopifnot(all(is.na(summary(result)$Minimum[, 3:5])))
result <- capture_warning(evolvabilityBetaMCMC2(G_mcmc, draws, post.dist = TRUE))
stopifnot(all(is.na(result$post.dist[, 3:5])),
          all(is.finite(result$post.dist[, 1:2])))

solve_g <- getFromNamespace("solve_g", "evolvability")
stopifnot(inherits(try(solve_g(matrix(1:6, 2, 3)), silent = TRUE), "try-error"))
