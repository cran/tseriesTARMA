#'  AR versus TARMA supLM robust test for nonlinearity
#'
#'  Implements a heteroskedasticity robust supremum Lagrange Multiplier test for a AR specification versus
#'  a TAR specification. Includes the classic (non robust) AR versus TAR test.
#'
#' @param x         A univariate time series.
#' @param pa        Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
#'                  The default is \code{0.25}
#' @param pb        Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
#'                  The default is \code{0.75}
#' @param ar.ord    Order of the AR part.
#' @param d         Delay parameter. Defaults to \code{1}.
#'
#' @details
#' Implements a heteroskedasticity robust asymptotic supremum Lagrange Multiplier test to test an AR specification versus a TAR specification.
#' This is an asymptotic test and the value of the test statistic has to be compared with the critical values tabulated 
#' in \insertCite{Gor21}{tseriesTARMA} or \insertCite{And03}{tseriesTARMA}.
#' Both the non-robust \code{supLM} and the robust \code{supLMh} statistics are returned.  
#'
#' @return
#'   An object of class \code{TARMAtest} with components:
#' \describe{
#'  \item{\code{statistic}}{A named vector with the values of the classic \code{supLM} and robust \code{supLMh} statistics.}
#'  \item{\code{parameter}}{A named vector: \code{threshold} is the value that maximises the Lagrange Multiplier values.}
#'  \item{\code{test.v}}{Matrix of values of the LM statistic for each threshold given in \code{thd.range}.}
#'  \item{\code{thd.range}}{Range of values of the threshold.}
#'  \item{\code{fit}}{The null model: AR fit over \code{x}.}
#'  \item{\code{sigma2}}{Estimated innovation variance from the AR fit.}
#'  \item{\code{data.name}}{A character string giving the name of the data.}
#'  \item{\code{prop}}{Proportion of values of the series that fall in the lower regime.}
#'  \item{\code{p.value}}{The p-value of the test. It is \code{NULL} for the asymptotic test.}
#'  \item{\code{method}}{A character string indicating the type of test performed.}
#'  \item{\code{d}}{The delay parameter.}
#'  \item{\code{pa}}{Lower threshold quantile.}
#'  \item{\code{dfree}}{Effective degrees of freedom. It is the number of tested parameters.}
#'}
#' @importFrom stats coef residuals as.ts
#' @export
#' @author Simone Giannerini, \email{simone.giannerini@@unibo.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Gor23}{tseriesTARMA}
#' * \insertRef{And03}{tseriesTARMA}
#'
#' @seealso \code{\link{TAR.test.B}} for the bootstrap version of the test. 
#' \code{\link{TARMA.test}} for the ARMA vs TARMA asymptotic version of the test, which includes also the AR vs TAR test, with different defaults.
#' \code{\link{TARMAGARCH.test}} for the robust version of the ARMA vs TARMA test with respect to 
#' GARCH innovations. 
#' \code{\link{TARMA.sim}} to simulate from a TARMA process. 
#'
#' @examples
#' set.seed(123)
#' ## a TAR(1,1) ---------------
#' x1   <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0.0,0.8), theta1=0, theta2=0, d=1, thd=0.2)
#' TAR.test(x1, ar.ord=1, d=1)
#'
#' ## a AR(1)    ----------------
#' x2   <- arima.sim(n=100, model=list(order=c(1,0,0), ar=0.5))
#' TAR.test(x2, ar.ord=1, d=1)
#'
#'
## ***************************************************************************

TAR.test <- function(x,pa=.25, pb=.75, ar.ord, d=1){

    # full Fortran version
    DNAME  <- deparse(substitute(x))
    n      <- length(x)
    k      <- max(ar.ord,d)
    neff   <- n-k
    a      <- ceiling((neff-1)*pa)
    b      <- floor((neff-1)*pb)
    nr     <- b-a+1
    test.v <- matrix(0,nr,2)
    storage.mode(test.v) <- 'double'
    # SUBROUTINE ARvsTARh(x,n,d,p,a,b,neff,testv,xth,coef,s2,eps)
    res    <- .Fortran('arvstarh',as.double(x),as.integer(n),as.integer(d),
    as.integer(ar.ord),as.integer(a),as.integer(b),as.integer(neff),
    test=test.v,xth=double(neff),coef=double(ar.ord+1),s2=double(1),eps=double(neff),PACKAGE='tseriesTARMA')
    test.v     <- res$test
    xth        <- res$xth
    names(res$coef) <- c('Intercept',paste('ar',1:ar.ord,sep=''))
    fit        <- list(coefficients=res$coef, sigma2=res$s2,residuals=res$eps)
    thd.range  <- sort(xth)[a:b]
    thd        <- thd.range[which.max(test.v[,1])]
    names(thd) <- 'threshold'
    names(d)   <- 'delay'
    test.stat  <- apply(test.v,FUN=max,MARGIN=2)
    names(test.stat) <- c('supLM','supLMh')
    colnames(test.v) <- c('supLM','supLMh')
    Il         <- (xth<= thd)
    dfree <- 1 + ar.ord
     METHOD <- paste('supLM test AR vs TAR. Null model: AR(',ar.ord,')',sep='')
    structure(list(test.v=test.v,thd.range=thd.range, fit=fit, sigma2=res$s2,
     parameter=c(thd), data.name=DNAME,prop=mean(Il), statistic=test.stat, p.value=NULL,
     method=METHOD,d=d,pa=pa, dfree=dfree),class=c('TARMAtest','htest'))
}

## ****************************************************************************

