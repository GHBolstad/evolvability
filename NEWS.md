# evolvability 2.0.2

## Bug fixes

* The check for symmetric G matrices never fired in `evolvabilityBeta()`,
  `evolvabilityMeans()`, `responsediffBeta()` and `conditionalG()` (it compared
  the upper triangle with itself). Asymmetric matrices now give an error in all
  functions that take a G matrix, including `meanStdG()`, `meanStdGMCMC()` and
  the MCMC functions (where every draw is checked), while rounding errors are
  tolerated.
* `macro_pred()` (and thereby `rate_gls(model = "recent_evol")` with
  `useLFO = TRUE`) gave wrong leave-one-out means when `y` contained tied
  values.
* `round_and_format()` now works on vectors when `sign_digits` is used.
* `Almer_SE()` now passes `control` on to every iteration and only warns about
  the maximum number of iterations when it has actually not converged.
* `Almer()` gives a clear error when a name in `A` does not match a random
  effect term.
* Singular `G` matrices are now detected with `rcond()` rather than by matching
  the text of an English error message.
* G matrices that cannot be inverted by `solve()` are now handled in the same
  way in all functions: a warning is given, the conditional evolvability,
  autonomy and integration are `NA`, and the other parameters are still
  calculated. This was not the case for `evolvabilityMeans()`, which silently
  returned `NaN`, `evolvabilityMeansMCMC()`, which failed with a cryptic error
  when a draw was singular, and `conditionalG()`, which gave a LAPACK error
  (it now returns a matrix of `NA`).
* A warning is now given if G is not positive definite (it has negative
  eigenvalues) in `evolvabilityBeta()`, `evolvabilityBetaMCMC()`,
  `evolvabilityBetaMCMC2()`, `evolvabilityMeans()`, `evolvabilityMeansMCMC()`,
  `conditionalG()` and `responsediffBeta()`. The results are still calculated.
  Some of the examples used matrices that were not positive definite and have
  been changed.
* G matrices with missing or non-finite values and selection gradients of the
  wrong dimension now give clear error messages.
* `evolvabilityMeansMCMC()` printed "Posterior means" for what are posterior
  medians.

## Other changes

* `evolvabilityMeans()` and `evolvabilityMeansMCMC()` compute the eigenvalues
  of `G` once and share their implementation (much faster for many draws).
* `means` must now have the same length as the dimension of `G` (or be `1`).
* `conditionalG()` validates `condition_on`.
* `GLS()` gives an informative error when neither `R` nor `L` is provided and
  documents that the first column of `X` must be the intercept.
* `rate_gls()` and `rate_gls_boot()` report progress with `message()` instead
  of `print()`, and obsolete commented-out code has been removed.
* The tests now use testthat, and tests were added for `GLS()`, `macro_pred()`,
  `conditionalG()`, `evolvabilityMeans()`, `round_and_format()` and
  `rate_gls()`.
* `DESCRIPTION` gained `URL`, `BugReports` and a minimal R version.
