#' @title multiMNbreaks: Structural Changes in Multivariate Regressions
#'
#' @description Estimation, testing and inference for multiple structural
#' changes at unknown dates in a system of equations, following Qu and Perron
#' (2007). Changes may occur in the regression coefficients, in the covariance
#' matrix of the errors, or in both, and linear restrictions may be imposed
#' within and across regimes. The package provides the quasi-maximum-likelihood
#' estimates of the break dates and coefficients, the sup-LR, WDmax and
#' sequential tests for the number of breaks, and simulated critical values
#' for any number of changing parameters.
#'
#' @details The computational core follows the GAUSS code distributed by the
#' authors: \code{\link{transf}} stacks the data, \code{\link{est}} estimates
#' the break dates and coefficients, \code{\link{suplr}}, \code{\link{supft}}
#' and \code{\link{seqlr}} compute the tests, and
#' \code{\link{critical_values}} supplies critical values.
#'
#' @references Qu, Z. and P. Perron (2007), "Estimating and Testing Structural Changes in Multivariate Regressions", \emph{Econometrica}, 75, 459-502.
#' @references Bai, J. and P. Perron (1998), "Estimating and Testing Linear Models with Multiple Structural Changes", \emph{Econometrica}, 66, 47-78.
#' @references Bai, J. and P. Perron (2003), "Computation and Analysis of Multiple Structural Change Models", \emph{Journal of Applied Econometrics}, 18, 1-22.
#'
#' @keywords internal
"_PACKAGE"
