#' @title Methods for TARMA tests
#'
#' @param x     A \code{TARMAtest} object.
#' @param \dots Further parameters passed to \code{\link{print.htest}}.
#' 
#' @details
#' Print method for \code{TARMAtest} objects. Prints the results using the method for class \code{htest} 
#' and adds critical values extracted from \code{\link{ACValues}} for the test for threshold nonlinearity
#' and \code{\link{supLMQur}} for the unit root test against the TARMA alternative.
#' Note that the bootstrap version of the tests also print the bootstrap p-value. 
#' @return
#'   No return value, called for side effects
#' 
#' @method print TARMAtest
#' @export
#' @author Simone Giannerini, \email{simone.giannerini@@uniud.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' 
#' @seealso \code{\link{print.htest}}

print.TARMAtest <- function(x,...){
  y <- x
  class(y) <- 'htest'
  print(y,...)
  if(is.null(x$dfree)){ # Unit Root tests
    cval <- round(CritValuesUR(x$fit.ARMA$nobs,x$pa,x$fit.ARMA$coef['ma1']),2)
    cat('delay: d = ',as.integer(x$d), ', lower threshold quantile: ',x$pa,'\n',sep='')
    cat('---------------      90%    95%    99%   99.9%   \n')
  }else{ # Classic tests
    cval <- CritValues(x$pa,x$dfree)
    cat('delay: d = ',as.integer(x$d),', degrees of freedom: ',as.integer(x$dfree),
        ', lower threshold quantile: ',x$pa,'\n',sep='')
    cat('---------------\n')
    cat('Level:               90%    95%    99%  \n')
  }
  cat('Critical Values: ',format(cval,width=6),'\n')
  cat("\n")
}
