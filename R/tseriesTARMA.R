#' tseriesTARMA: Analysis of Nonlinear Time Series through Threshold Autoregressive Moving Average Models (TARMA) models
#'
#' It provides advanced functions for:
#' * TARMA model fitting and forecasting:
#'    - Least Squares fitting of a full subset TARMA model, including robust LS fitting.
#'    - Maximum Likelihood fitting of a subset TARMA model with common MA parts and possible covariates.
#' * TARMA testing for threshold type nonlinearity:
#'    - Tests for AR vs TAR (asymptotic, bootstrap, wild bootstrap)
#'    - Tests for ARMA vs TARMA with both i.i.d. errors and GARCH errors.
#' * Unit-root testing against a stationary TARMA model
#'
#' @author Simone Giannerini, \email{simone.giannerini@@unibo.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Gia23}{tseriesTARMA}
#' * \insertRef{Gia22}{tseriesTARMA}
#' * \insertRef{Gia21}{tseriesTARMA}
#' * \insertRef{Gor21}{tseriesTARMA}
#' * \insertRef{Gor23}{tseriesTARMA}
#' * \insertRef{Cha19}{tseriesTARMA}
#' * \insertRef{Cha24}{tseriesTARMA}
#' 
#' @docType package
#' @name tseriesTARMA
#' @useDynLib tseriesTARMA
#' @importFrom Rdpack reprompt
#' @import mathjaxr
#' @keywords internal
"_PACKAGE"
NULL
#> NULL
