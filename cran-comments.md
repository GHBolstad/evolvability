## Resubmission

This is a resubmission of a package that was archived on 2026-11-12 because
issues were not corrected in time. The issues were caused by the update of
lme4 (the way weights and other arguments are evaluated in `lme4::lFormula()`),
which broke `Almer()`. In this version (2.0.2), `Almer()` passes evaluated
arguments to `lFormula()`, and regression tests for weights, offsets and
formula environments have been added. See `NEWS.md` for the other changes.

## Test environments

* Local: Ubuntu 22.04, R 4.4.1, with lme4 2.1.0 (the current CRAN version)

## R CMD check results

`R CMD check --as-cran` gives 0 errors | 0 warnings | 2 notes:

* "New submission" and "Package was archived on CRAN": see above.
* "unable to verify current time": a quirk of the local check environment.

## Downstream dependencies

There are no reverse dependencies.
