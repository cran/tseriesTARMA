#' Tabulated Critical Values for the Unit Root IMA vs TARMA test 
#'
#' @description \loadmathjax
#' The data provides asymptotic null critical 
#' values for the unit root supLM test described in \insertCite{Cha24}{tseriesTARMA}. 
#' They are used with the following tests: 
#' - \code{\link{TARMAur.test}}
#' - \code{\link{TARMAur.test.B}}
#' 
#' Provided \code{pb = 1- pa}.
#' 
#' @format ## `supLMQur`
#' A 4-dimensional array that contains 4 critical values (at level 10%, 5%, 1%, 0.1%)
#' for each combination of
#' \describe{
#'   \item{\code{pa}}{Lower bound for the threshold range. From 0.01 to 0.4}
#'   \item{\code{th}}{MA(1) parameter.}
#'   \item{\code{n}}{Sample size.}
#' }
#' @source \insertRef{Cha24}{tseriesTARMA}
"supLMQur"