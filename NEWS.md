# mvbreak 0.0.1.9000 (development)

First usable version. The computational core is a port of the authors' GAUSS
code, validated against it on the two data sets it ships with; everything
user-facing is new.

## Model functions

* `breaks_varp()`, `breaks_mvreg()`, `breaks_sur()`, `breaks_panel()` and the
  general `breaks_system()` fit a system with `m` breaks. Arguments `change`
  (`"coef"`, `"cov"`, `"both"`), `fixed` and `fixed_equations` replace the
  flags and hand-built restriction matrices of the original code.
* `breaks_ordered()` fits the locally ordered breaks model (Section 5 of the
  paper) for a multivariate regression, a SUR or a VAR.

## Tests, selection and intervals

* `test_breaks()` computes the sup-LR (or HAC sup-F) tests for 1 to
  `max_breaks` breaks, UDmax, WDmax and the sequential tests, with critical
  values; `select_breaks()` chooses the number of breaks sequentially.
* `critical_values()` reads a simulated table covering up to 60 changing
  parameters, five trimmings and four levels, and simulates anything outside it.
* `confint_breaks()` gives the break-date intervals of Theorem 3, with
  options for serial correlation, prewhitening and regime-specific
  regressors, and the simulated intervals of Theorem 4 for ordered breaks.

## Methods and simulation

* `print`, `summary`, `coef`, `vcov`, `residuals`, `fitted`, `logLik`, `nobs`,
  `confint` and `plot` for fitted models; `summary` for test objects.
* `simulate_varp()`, `simulate_mvreg()`, `simulate_sur()`, `simulate_panel()`
  and `simulate_ordered()` generate data with known breaks.

## Differences from the original GAUSS code

* The sequential test trims each segment relative to its own length (Section
  6.3 of the paper) and imposes a floor on the regime length so that every
  covariance matrix is nonsingular; it returns `NA`, not 0, when no segment
  can hold another break.
* The prewhitened long-run variance in the intervals is computed on the
  distinct elements of the residual products, which removes a singularity;
  the quantile search is exact, with an adaptive bracket.
* The trimming rounds `trim * T` with R's `round()`; `h` can be passed to
  reproduce the original rounding at half-integers. Intervals are clipped to
  the sample.
* The intervals for ordered breaks use the residuals of the fitted model;
  the original code used those of the model without breaks.
* Critical values for the restricted model use the fitted number of breaks.
* Dense `nT x nT` matrices are replaced by sums over observations; the
  ported internals are not exported.
