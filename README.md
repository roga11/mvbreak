# mvbreak

Structural changes at unknown dates in systems of equations, following
Qu and Perron (2007, *Econometrica*): estimation, tests for the number of
breaks, and confidence intervals for the break dates, with changes in the
regression coefficients, in the covariance matrix of the errors, or in both,
and with coefficients that can be held fixed across regimes. The package
also implements the locally ordered breaks model of Section 5 of the paper.

The package is under development and its name is provisional.

## Installation

```r
# install.packages("pak")
pak::pak("roga11/mvbreak")
# or: remotes::install_github("roga11/mvbreak")
```

## A first example

The two-equation system of the original GAUSS example: each equation has an
intercept and its own regressor, and two breaks are allowed in every
coefficient.

```r
library(mvbreak)
y <- as.matrix(data1[, 1:2]); z <- as.matrix(data1[, 3:5])

fit <- breaks_sur(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), m = 2, trim = 0.2)
fit$dates                 # 56, 100
summary(fit)              # coefficients by regime with standard errors

test_breaks(fit)          # sup-LR, UDmax, WDmax and sequential tests with critical values
confint_breaks(fit)       # intervals for the break dates
plot(fit, ci = confint_breaks(fit))

select_breaks(breaks_sur(y, x = list(z[, 2, drop = FALSE], z[, 3, drop = FALSE]), m = 3, trim = 0.2))
```

## Functions

| Purpose | Functions |
|---|---|
| Fit a model with `m` breaks | `breaks_varp()`, `breaks_mvreg()`, `breaks_sur()`, `breaks_panel()`, `breaks_system()` |
| Locally ordered breaks | `breaks_ordered()` |
| Tests and the number of breaks | `test_breaks()`, `select_breaks()`, `critical_values()` |
| Confidence intervals for the dates | `confint_breaks()`, `confint()` |
| Methods on a fit | `print`, `summary`, `coef`, `vcov`, `residuals`, `fitted`, `logLik`, `nobs`, `plot` |
| Simulation | `simulate_varp()`, `simulate_mvreg()`, `simulate_sur()`, `simulate_panel()`, `simulate_ordered()` |

Every model function takes `m` (number of breaks), `trim` (minimum regime
length as a fraction of the sample), `change` (`"coef"`, `"cov"` or `"both"`)
and `fixed` (coefficients held constant across regimes, by name). The
vignette walks through all of them.

## Validation

The package was written from the authors' GAUSS code and checked against it:
on the two data sets distributed with that code, over 324 specifications
(model type, what changes, serial correlation, number of breaks, trimming),
every estimate, test statistic and confidence interval reproduces the
reference or differs for a documented reason. Critical values are read from a
table simulated for up to 60 changing parameters per break, checked against
Bai and Perron (1998, 2003), and simulated beyond it.

## References

Bai, J. and P. Perron (1998), "Estimating and Testing Linear Models with
Multiple Structural Changes", *Econometrica*, 66, 47–78.

Bai, J. and P. Perron (2003), "Critical Values for Multiple Structural Change
Tests", *Econometrics Journal*, 6, 72–78.

Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in
Multivariate Regressions", *Econometrica*, 75, 459–502.
