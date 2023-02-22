#'  Unit root supLM test for an integrated MA versus a stationary TARMA process
#'
#'  Implements a supremum Lagrange Multiplier unit root test for the null hypothesis of a integrated MA  process versus
#'  a stationary TARMA process.
#'
#' @param x         A univariate vector or time series.
#' @param B         Integer. Number of bootstrap resamples. Defaults to 1000.
#' @param pa        Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
#'                  The default is \code{0.25}
#' @param pb        Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
#'                  The default is \code{0.75}
#' @param thd.range Vector of optional user defined threshold range. If missing then \code{pa} and \code{pb} are used.
#' @param method    Fitting method to be passed to \code{arima}.
#' @param btype     Bootstrap type, can be one of \code{'iid','wb.h','wb.r','wb.n'}, see Details.
#' @param \dots     Additional arguments to be passed to \code{arima}.
#'
#' @details
#' Implements the bootstrap version of \code{\link{TARMAur.test}} the supremum Lagrange Multiplier test 
#' to test an integrate MA(1) specification versus a stationary TARMA(1,1) specification.
#' The option \code{btype} specifies the type of bootstrap as follows:
#'  \describe{
#'  \item{\code{wb.r}}{Residual wild bootstrap with Rademacher auxiliary distribution. See \insertCite{Gia22}{tseriesTARMA}.}
#'  \item{\code{wb.n}}{Residual wild bootstrap with Normal auxiliary distribution. See \insertCite{Gia22}{tseriesTARMA}.}
#'  \item{\code{iid}}{Residual iid bootstrap. See \insertCite{Gor21b}{tseriesTARMA}.}
#' }
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
#'  \item{\code{p.value}}{The bootstrap p-value of the test.}
#'  \item{\code{method}}{A character string indicating the type of test performed.}
#'  \item{\code{d}}{The delay parameter.}
#'  \item{\code{pa}}{Lower threshold quantile.}
#'  \item{\code{Tb}}{The bootstrap null distribution.}
#'}
#' @importFrom stats coef residuals
#' @export
#' @author Simone Giannerini, \email{simone.giannerini@@unibo.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Cha24}{tseriesTARMA}
#'
#' @seealso \code{\link{TARMAur.test}} for the asymptotic version of the test. \code{\link{print.TARMAtest}} for the print method.
#'
#' @examples
#' ## a TARMA(1,1,1,1) 
#' set.seed(123)
#' x1    <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0.0,0.8), theta1=0.5, theta2=0.5, d=1, thd=0.2)
#' TARMAur.test.B(x1, B=100) # B=100 for speedup
#'
#'
#' ## a IMA(1,1)
#' x2   <- arima.sim(n=100, model=list(order = c(0,1,1),ma=0.6))
#' TARMAur.test.B(x2, B=100) # B=100 for speedup
#'
#'
## ***************************************************************************

TARMAur.test.B <- function(x, B=1000, pa=.25, pb=.75, thd.range, method='ML',
                              btype = c('wb.r','wb.n','iid'), ...){
  DNAME  <- deparse(substitute(x))
  if(!is.ts(x)) x <- ts(x)
  btype  <- match.arg(btype)
  n      <- length(x)
  fit    <- stats::arima(x, order=c(0,1,1), method=method,xreg=data.frame(intercept=1:n))
    if(fit$code!=0) stop('ARMA fit does not converge')
  # elements from ARMA fit that we need**************************************
  ma.1      <- -coef(fit)['ma1']
  epst      <- residuals(fit)
  resi      <- scale(epst,center=TRUE,scale=FALSE) # centered residuals for the bootstrap
  s2        <- fit$sigma2                          # taken as estimate of sigma^2
  # *****************************************************************************
  nstart <- 1
  np     <- n + nstart # transient
  if(btype == 'wb.r'){
    eta <- matrix(sample(c(-1,1),replace=TRUE, size=np*B),np,B)
    resi.b    <- c(rnorm(nstart),resi)*eta
    smeth     <- 'wild bootstrap (Rademacher)'
  }else
    if(btype == 'wb.n'){
      eta <- matrix(rnorm(np*B),np,B)
      resi.b    <- c(rnorm(nstart),resi)*eta
      smeth     <- 'wild bootstrap (Gaussian)'
    }else
      if(btype == 'iid'){
        resi.b    <- matrix(sample(resi,size=np*B,replace=TRUE),np,B)
        smeth     <- 'i.i.d. bootstrap'
    }
  # ******************************************************************************
  test.a <- TARMAur.test(x,pa=pa, pb=pb, method=method,...)
  
  xt     <- x[2:n];     # X_t
  xth    <- x[1:(n-1)]  # X_(t-1)
  neff   <- length(xt)  # effective sample size
  x.b    <- apply(rbind(rep(x[1],B),resi.b[-1,]-ma.1*resi.b[-n,]),MARGIN=2,FUN=cumsum)
  x.b    <- x.b[-(1:nstart),]
  Tb     <- double(B)
  for(k in 1:B){
    Tb[k]  <- TARMAur.test(x.b[,k], pa=pa, pb=pb, method=method, ...)$statistic
  }
  test.stat <- test.a$statistic
  names(test.stat) <- 'supLM'
  p.value   <- mean(Tb>test.stat,na.rm=TRUE)
  METHOD    <- paste(smeth,' supLM unit root test IMA vs TARMA',sep='')
  structure(list(statistic=test.a$statistic, parameter=test.a$parameter, 
    test.v=test.a$test.v, thd.range=test.a$thd.range, fit.ARMA=test.a$fit.ARMA, 
    sigma2=test.a$sigma2, data.name=DNAME, p.value=p.value, method=METHOD,
    d=1, pa=pa, Tb=Tb),class=c('TARMAtest','htest'))
}

### ***************************************************************************