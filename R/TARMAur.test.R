#'  Unit root supLM test for an integrated MA versus a stationary TARMA process
#'
#'  Implements a supremum Lagrange Multiplier unit root test for the null hypiythesis of a integrated MA  process versus
#'  a stationary TARMA process.
#'
#' @param x         A univariate vector or time series.
#' @param pa        Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
#'                  The default is \code{0.25}
#' @param pb        Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
#'                  The default is \code{0.75}
#' @param thd.range Vector of optional user defined threshold range. If missing then \code{pa} and \code{pb} are used.
#' @param method    Fitting method to be passed to \code{arima}.
#' @param \dots     Additional arguments to be passed to \code{arima}.
#'
#' @details
#' Implements an asymptotic supremum Lagrange Multiplier test to test an integrate MA(1) specification versus a 
#' stationary TARMA(1,1) specification. This is an asymptotic test and the value of the test statistic has to be compared with the critical 
#' values tabulated in \insertCite{Cha24}{tseriesTARMA} and available in \code{\link{supLMQur}}. 
#' The relevant critical values are automatically shown upon printing the test, see \code{\link{print.TARMAtest}}.
#'
#' @return
#'   An object of class \code{TARMAtest} with components:
#' \describe{
#'  \item{\code{statistic}}{The value of the supLM statistic.}
#'  \item{\code{parameter}}{A named vector: \code{threshold} is the value that maximises the Lagrange Multiplier values.}
#'  \item{\code{test.v}}{Vector of values of the LM statistic for each threshold given in \code{thd.range}.}
#'  \item{\code{thd.range}}{Range of values of the threshold.}
#'  \item{\code{fit.ARMA}}{The null model: IMA(1) fit over \code{x}.}
#'  \item{\code{sigma2}}{Estimated innovation variance from the IMA fit.}
#'  \item{\code{data.name}}{A character string giving the name of the data.}
#'  \item{\code{p.value}}{The p-value of the test. It is \code{NULL} for the asymptotic test.}
#'  \item{\code{method}}{A character string indicating the type of test performed.}
#'  \item{\code{d}}{The delay parameter.}
#'  \item{\code{pa}}{Lower threshold quantile.}
#'}
#' @importFrom stats coef residuals
#' @export
#' @author Simone Giannerini, \email{simone.giannerini@@uniud.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Cha24}{tseriesTARMA}
#'
#' @seealso \code{\link{TARMAur.test.B}} for the bootstrap version of the test. 
#' \code{\link{print.TARMAtest}} for the print method.
#'
#' @examples
#' ## a TARMA(1,1,1,1) 
#' set.seed(123)
#' x1    <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0.0,0.8), theta1=0.5, theta2=0.5, d=1, thd=0.2)
#' TARMAur.test(x1)
#'
#'
#' ## a IMA(1,1)
#' x2   <- arima.sim(n=100, model=list(order = c(0,1,1),ma=0.6))
#' TARMAur.test(x2)
#'
#'
## ***************************************************************************

TARMAur.test <- function(x,pa=.25, pb=.75, thd.range, method='ML',...){
    
    DNAME  <- deparse(substitute(x))
    if(!is.ts(x)) x <- ts(x)
    n      <- length(x)
    fit    <- stats::arima(x, order=c(0,1,1), method=method,xreg=data.frame(intercept=1:n))
    thd.v  <- sort(x)
    a      <- ceiling((n-1)*pa)
    b      <- floor((n-1)*pb)
    ma     <- -coef(fit)['ma1']      # notice the change in sign
    epst   <- residuals(fit) # taken as estimates of eps_{t},
    s2     <- fit$sigma2
    if(missing(thd.range)) thd.range <- thd.v[a:b]
    x1     <- x[-n]              # X_{t-1}
    I1     <- rep(1,n-1)         # intercept
    nr     <- length(thd.range)  # number of points for threshold estimation
    
    # SUBROUTINE IMAvsTARMA(x,eps,n,trange,nr,s2,ma,testv)
    res    <- .Fortran('imavstarma',as.double(x),as.double(epst),as.integer(n),
                       as.double(thd.range),as.integer(nr),as.double(s2),as.double(ma),test=double(nr),PACKAGE='tseriesTARMA')
    test.v           <- res$test
    thd              <- thd.range[which.max(test.v)]
    names(thd)       <- 'threshold'
    test.stat        <- max(test.v)
    names(test.stat) <- 'supLM'
    p.value          <- NULL
    METHOD <- paste('supLM unit root test IMA vs TARMA',sep='')
    structure(list(statistic=test.stat, parameter=thd, 
      test.v=test.v, thd.range=thd.range, fit.ARMA=fit, 
      sigma2=s2, data.name=DNAME, p.value=p.value, method=METHOD,
      d=1, pa=pa),class=c('TARMAtest','htest'))
  }
  
  ## ***************************************************************************  
