#'  ARMA versus TARMA supLM test for nonlinearity
#'
#'  Implements a supremum Lagrange Multiplier test for a ARMA specification versus
#'  a TARMA specification. Includes the AR versus TAR test.
#'
#' @param x         A univariate time series.
#' @param pa        Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
#'                  The default is \code{0.25}
#' @param pb        Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
#'                  The default is \code{0.75}
#' @param ar.ord    Order of the AR part.
#' @param ma.ord    Order of the MA part.
#' @param ma.fixed  Logical. Only applies to testing ARMA vs TARMA. If \code{TRUE} computes the test where only the AR parameters are tested, see \insertCite{Gor21}{tseriesTARMA} for details.
#' @param d         Delay parameter. Defaults to \code{1}.
#' @param thd.range Vector of optional user defined threshold range. If missing then \code{pa} and \code{pb} are used.
#' @param method    Fitting method to be passed to \code{arima}.
#' @param \dots     Additional arguments to be passed to \code{arima}.
#'
#' @details
#' Implements an asymptotic supremum Lagrange Multiplier test to test an ARMA specification versus a TARMA specification.
#' If \code{ma.fixed=TRUE} (the default), the AR parameters are tested whereas the MA parameters are fixed. If \code{ma.fixed=FALSE} both the AR and the MA parameters are tested.
#'  This is an asymptotic test and the value of the test statistic has to be compared with the critical values tabulated in \insertCite{Gor21}{tseriesTARMA} and \insertCite{And03}{tseriesTARMA}.
#' If \code{ma.ord=0} then the AR versus TAR test is used. Note that when \code{method='CSS'}, this is equivalent to \code{TAR.test}, which uses least squares.
#'
#' @return
#'   An object of class \code{TARMAtest} with components:
#' \describe{
#'  \item{\code{statistic}}{The value of the supLM statistic.}
#'  \item{\code{parameter}}{A named vector: \code{threshold} is the value that maximises the Lagrange Multiplier values.}
#'  \item{\code{test.v}}{Vector of values of the LM statistic for each threshold given in \code{thd.range}.}
#'  \item{\code{thd.range}}{Range of values of the threshold.}
#'  \item{\code{fit.ARMA}}{The null model: ARMA fit over \code{x}.}
#'  \item{\code{sigma2}}{Estimated innovation variance from the ARMA fit.}
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
#' @seealso \code{\link{TAR.test.B}} for the bootstrap version of the test. \code{\link{TARMAGARCH.test}}
#' for the robust version of the test with respect to GARCH innovations. \code{\link{TARMA.sim}} to simulate from a TARMA process.
#'
#' @examples
#' ## a TARMA(1,1,1,1) where the threshold effect is on the AR parameters
#' set.seed(123)
#' x1    <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0.0,0.8), theta1=0.5, theta2=0.5, d=1, thd=0.2)
#' TARMA.test(x1, ar.ord=1, ma.ord=1, d=1)
#' TARMA.test(x1, ar.ord=1, ma.ord=1, d=1, ma.fixed=FALSE) # full TARMA test
#'
#' ## a TARMA(1,1,1,1) where the threshold effect is on the MA parameters
#' set.seed(212)
#' x2    <- TARMA.sim(n=100, phi1=c(0.5,0.2), phi2=c(0.5,0.2), theta1=0.6, theta2=-0.6, d=1, thd=0.2)
#' TARMA.test(x2, ar.ord=1, ma.ord=1, d=1)
#' TARMA.test(x2, ar.ord=1, ma.ord=1, d=1, ma.fixed=FALSE) # full TARMA test
#'
#' ## a ARMA(1,1)
#' x3   <- arima.sim(n=100, model=list(order = c(1,0,1),ar=0.5, ma=0.5))
#' TARMA.test(x3, ar.ord=1, ma.ord=1, d=1)
#'
#' ## a TAR(1,1)
#' x4   <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0.0,0.8), theta1=0, theta2=0, d=1, thd=0.2)
#' TARMA.test(x4, ar.ord=1, ma.ord=0, d=1)
#'
#' ## a AR(1)
#' x5   <- arima.sim(n=100, model=list(order = c(1,0,0),ar=0.5))
#' TARMA.test(x5, ar.ord=1, ma.ord=0, d=1)
#'
#'
## ***************************************************************************'

TARMA.test <- function(x,pa=.25, pb=.75, ar.ord, ma.ord, ma.fixed=TRUE, d, thd.range, method='CSS-ML',...){

    DNAME  <- deparse(substitute(x))
    n      <- length(x)
    # thd range **********************************
    if(missing(thd.range)){
        thd.v     <- sort(x[-1])
        a         <- ceiling((n-1)*pa)
        b         <- floor((n-1)*pb)
        thd.range <- thd.v[a:b]
    }
    # ********************************************
    nr   <- length(thd.range)   # number of points for threshold estimation
    fit    <- stats::arima(x, order=c(ar.ord,0,ma.ord), method=method) #fit ARMA(p,q)
    if(fit$code!=0) stop('ARMA fit does not converge')
    if(ma.ord>0){
        ma     <- -coef(fit)[-c(1:ar.ord,ar.ord+ma.ord+1)]  # notice the change in sign
        if(any(Mod(polyroot(c(1, -ma))) < 1.10)){
            fit    <- stats::arima(x, order=c(ar.ord,0,ma.ord), method='CSS') #fit ARMA(p,q) with CSS
            if(fit$code!=0) stop('ARMA fit does not converge')
            ma     <- -coef(fit)[-c(1:ar.ord,ar.ord+ma.ord+1)]  # notice the change in sign
            if(any(Mod(polyroot(c(1, -ma))) < 1.04)) stop('non invertible MA part')
        }
     }
    # elements from ARMA fit that we need*********************************************
    epst   <- residuals(fit)                            # taken as estimates of eps_{t}
    s2     <- fit$sigma2                                # taken as estimate of sigma2
    # *********************************************************************************

    test.v <- double(nr)
    if(ma.ord==0){
        res    <- .Fortran('arvstar',as.double(x),as.double(epst),as.integer(n),as.integer(d),
        as.double(thd.range),as.integer(nr),as.double(s2),as.integer(ar.ord),test=test.v, PACKAGE='tseriesTARMA')
    }else
    if(ma.ord>0){
        if(ma.fixed){
        # SUBROUTINE ARMAvsTARMA(x,eps,n,d,trange,nr,s2,p,ma,q,testv)
        res    <- .Fortran('armavstarma',as.double(x),as.double(epst),as.integer(n),as.integer(d),
        as.double(thd.range),as.integer(nr),as.double(s2),as.integer(ar.ord),as.double(ma),
        as.integer(ma.ord),test=test.v,PACKAGE='tseriesTARMA')
        }else{
        # SUBROUTINE ARMAvsTARMAg(x,eps,n,d,trange,nr,s2,p,ma,q,testv)
        res    <- .Fortran('armavstarmag',as.double(x),as.double(epst),as.integer(n),as.integer(d),
        as.double(thd.range),as.integer(nr),as.double(s2),as.integer(ar.ord),as.double(ma),
        as.integer(ma.ord),test=test.v,PACKAGE='tseriesTARMA')
        }
    }
    test.v     <- res$test
    thd        <- thd.range[which.max(test.v)]
    names(thd) <- 'threshold'
    names(d)   <- 'delay'
    test.stat  <- max(test.v)
    names(test.stat) <- 'supLM'
    x          <- as.ts(x) # this makes ts.intersect work with zoo objects
    xreg       <- ts.intersect(x,lag(x,-d))
    xth        <- xreg[,2]
    Il         <- (xth<= thd)
    if(ma.ord==0){
        METHOD <- paste('supLM test AR vs TAR. Null model: AR(',ar.ord,')',sep='')
    }else{
        fstring <- ifelse(ma.fixed,'AR','both AR and MA')
        METHOD  <- paste('supLM test ARMA vs TARMA on ',fstring,' parameters. Null model: ARMA(',ar.ord,',',ma.ord,')',sep='')
    }
    dfree <- 1 + ar.ord + ma.ord*(!ma.fixed) # effective degrees of freedom (number of tested parameters)
    structure(list(test.v=test.v,
     thd.range=thd.range, fit.ARMA=fit, sigma2=s2,
     parameter=c(thd), data.name=DNAME,prop=mean(Il), statistic=test.stat, p.value=NULL,
     method=METHOD,d=d,pa=pa, dfree=dfree),class=c('TARMAtest','htest'))
}

## ****************************************************************************
