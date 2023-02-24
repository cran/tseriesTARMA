#'  ARMA GARCH versus TARMA GARCH supLM test for nonlinearity
#'
#'  Implements a supremum Lagrange Multiplier test for a ARMA-GARCH specification versus
#'  a TARMA-GARCH specification. Both the AR and MA parameters are tested
#'
#' @param x         A univariate time series.
#' @param pa        Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
#'                  The default is \code{0.25}
#' @param pb        Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
#'                  The default is \code{0.75}
#' @param ar.ord    Order of the AR part.
#' @param ma.ord    Order of the MA part.
#' @param arch.ord  Order of the ARCH part.
#' @param garch.ord Order of the GARCH part.
#' @param d         Delay parameter. Defaults to \code{1}.
#' @param thd.range Vector of optional user defined threshold range. If missing then \code{pa} and \code{pb} are used.
#' @param method    Fitting method to be passed to \code{arima}.
#' @param \dots     Additional arguments to be passed to \code{arima}.
#'
#' @details
#' Implements an asymptotic supremum Lagrange Multiplier test to test an ARMA-GARCH specification versus
#' a TARMA-GARCH specification. In other words, the test is robust with respect to heteroskedasticity.
#' Both the AR parameters and the MA parameters are tested. This is an asymptotic test and the value of
#' the test statistic has to be compared with the critical values tabulated in \insertCite{Ang22}{tseriesTARMA}
#'  or \insertCite{And03}{tseriesTARMA}.
#'
#' @return
#'   A list of class \code{htest} with components:
#' \describe{
#'  \item{\code{statistic}}{The value of the supLM statistic.}
#'  \item{\code{parameter}}{A named vector: \code{threshold} is the value that maximises the Lagrange Multiplier values.}
#'  \item{\code{test.v}}{Vector of values of the LM statistic for each threshold given in \code{thd.range}.}
#'  \item{\code{thd.range}}{Range of values of the threshold.}
#'  \item{\code{fit.ARMA}}{The null model: ARMA fit over \code{x}.}
#'  \item{\code{fit.GARCH}}{The null model: GARCH fit over the residuals of the ARMA fit.}
#'  \item{\code{sigma2}}{Estimated innovation variance from the ARMA fit.}
#'  \item{\code{data.name}}{A character string giving the name of the data.}
#'  \item{\code{prop}}{Proportion of values of the series that fall in the lower regime.}
#'  \item{\code{p.value}}{The p-value of the test. It is \code{NULL} for the asymptotic test.}
#'  \item{\code{method}}{A character string indicating the type of test performed.}
#'  \item{\code{d}}{The delay parameter.}
#'  \item{\code{pa}}{Lower threshold quantile.}
#'  \item{\code{dfree}}{Effective degrees of freedom. It is the number of tested parameters.}
#'}
#' @importFrom stats coef residuals
#' @importFrom rugarch ugarchspec ugarchfit
#' @export
#' @author Simone Giannerini, \email{simone.giannerini@@unibo.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Ang22}{tseriesTARMA}
#' * \insertRef{Gor21}{tseriesTARMA}
#' * \insertRef{And03}{tseriesTARMA}
#'
#' @seealso \code{\link{TARMA.test}} and \code{\link{TAR.test.B}} for the asymptotic and bootstrap test without the GARCH component. \code{\link{TARMA.sim}} to simulate from a TARMA process. \code{\link{TARMA.fit}} and \code{\link{TARMA.fit2}} for TARMA modelling.
#'
#' @examples
#' ## Function to simulate from a ARMA-GARCH process
#'
#' arma11.garch11 <- function(n, ph, th, a, b, a0=1, rand.gen = rnorm, innov = rand.gen(n, ...),
#' n.start = 500, start.innov = rand.gen(n.start, ...),...){
#'
#'   #  Simulates a ARMA(1,1)-GARCH(1,1) process
#'   #  with parameeters ph, th, a, b, a0.
#'   #         x[t] <- ph*x[t-1] + th*eps[t-1] + eps[t]
#'   #       eps[t] <- e[t]*sqrt(v[t])
#'   #         v[t] <- a0 + a*eps[t-1]^2 + b*v[t-1];
#'   # ph  : AR
#'   # th  : MA
#'   # a   : ARCH
#'   # b   : GARCH
#'
#'   # checks
#'   if(abs(a+b)>=1)   stop("model is not stationary")
#'   if(b/(1-a)>=1) stop("model has infinite fourth moments")
#'
#'   if (!missing(start.innov) && length(start.innov) < n.start)
#'     stop(gettextf("'start.innov' is too short: need %d points", n.start), domain = NA)
#'   e <- c(start.innov[1L:n.start], innov[1L:n])
#'   ntot <- length(e)
#'   x <- v <- eps <- double(ntot)
#'   v[1]   <- a0/(1.0-a-b);
#'   eps[1] <- e[1]*sqrt(v[1])
#'   x[1]   <- eps[1];
#'   for(i in 2:ntot){
#'     v[i]   <- a0 + a*eps[i-1]^2 + b*v[i-1];
#'     eps[i] <- e[i]*sqrt(v[i])
#'     x[i]   <- ph*x[i-1] + th*eps[i-1] + eps[i]
#'   }
#'   if (n.start > 0)  x <- x[-(1L:n.start)]
#'   return(ts(x));
#' }
#'
#' ## **************************************************************************
#' ## Comparison between the robust and the non-robust test in presence of GARCH errors
#' ## Simulates from the ARMA(1,1)-GARCH(1,1)
#'
#' set.seed(12)
#' x1 <- arma11.garch11(n=100, ph=0.9, th=0.5, a=0.85, b=0.1, a0=1,n.start=500)
#' TARMAGARCH.test(x1, ar.ord=1, ma.ord=1, arch.ord=1, garch.ord=1, d=1)
#' TARMA.test(x1, ar.ord=1, ma.ord=1, d=1, ma.fixed=FALSE)
#'
#' ## a TARMA(1,1,1,1) where the threshold effect is on the AR parameters
#' set.seed(123)
#' x2  <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0.0,0.8), theta1=0.5, theta2=0.5, d=1, thd=0.2)
#' TARMAGARCH.test(x2, ar.ord=1, ma.ord=1, d=1)
#'
#'
## ***************************************************************************

TARMAGARCH.test <- function(x,pa=.25, pb=.75, ar.ord=1, ma.ord=1, arch.ord=1, garch.ord=1, d=1, thd.range,  method='CSS', ...){
  DNAME  <- deparse(substitute(x))
  if(!is.ts(x)) x <- ts(x)
  n     <- length(x)
  shift <- mean(x)
  x     <- scale(x,scale=FALSE)
  # thd range **********************************
  if(missing(thd.range)){
    thd.v     <- sort(x[-1])
    a         <- ceiling((n-1)*pa)
    b         <- floor((n-1)*pb)
    thd.range <- thd.v[a:b]
  }
  nr   <- length(thd.range)   # number of points for threshold estimation
  # ********************************************

  fit    <- stats::arima(x, order=c(ar.ord,0,ma.ord),method=method,...) #fit ARMA(p,q)
  if(fit$code!=0) stop('ARMA fit does not converge')
  ma     <- -coef(fit)[-c(1:ar.ord,ar.ord+ma.ord+1)]  # notice the change in sign
  # if(any(Mod(polyroot(c(1, -ma))) < 1.05)){
  #   fit    <- stats::arima(x, order=c(ar.ord,0,ma.ord), method='CSS') #fit ARMA(p,q) with CSS
  #   if(fit$code!=0) stop('ARMA fit does not converge')
  #   ma     <- -coef(fit)[-c(1:ar.ord,ar.ord+ma.ord+1)]  # notice the change in sign
  #   if(any(Mod(polyroot(c(1, -ma))) < 1.04)) stop('non invertible MA part')
  # }

  # elements from ARMA fit that we need*********************************************
  epst     <- residuals(fit)                            # taken as estimates of eps_{t}
  s2       <- fit$sigma2                                # taken as estimate of sigma^2
  # *********************************************************************************
  # rugarch
  spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(arch.ord, garch.ord)),
                      mean.model = list(armaOrder=c(0,0),include.mean = FALSE))
  fit.g <- ugarchfit(spec=spec1, data = epst)
  ht    <- ts(fit.g@fit$var) # conditional variance
  aa    <- fit.g@fit$coef[2:(arch.ord+1)] # ARCH pars (without a0)
  bb    <- fit.g@fit$coef[(arch.ord+2):(arch.ord+garch.ord+1)] # GARCH pars
  if(fit.g@fit$convergence != 0){warning('GARCH fit did not reach convergence.')}
  ## call to Fortran code
  test.v <- double(nr)
  #SUBROUTINE TARMAGARCH(x,eps,h,n,d,trange,nr,p,ma,q,aa,m,bb,s,testv)
  res    <- .Fortran('tarmagarch',as.double(x),as.double(epst),as.double(ht),
                     as.integer(n),as.integer(d), as.double(thd.range),as.integer(nr),as.integer(ar.ord),
                     as.double(ma), as.integer(ma.ord),as.double(aa), as.integer(arch.ord),
                     as.double(bb), as.integer(garch.ord),test=test.v,PACKAGE='tseriesTARMA')

  test.v     <- res$test
  thd        <- thd.range[which.max(test.v)] + shift
  names(thd) <- 'threshold'
  names(d)   <- 'delay'
  test.stat  <- max(test.v)
  names(test.stat) <- 'supLM'
  xreg       <- ts.intersect(x,lag(x,-d))
  xth        <- xreg[,2]
  Il         <- (xth<= thd)
  dfree      <- 1 + ar.ord + ma.ord
  if(ma.ord==0){
    METHOD <- paste('supLM test AR-GARCH vs TAR-GARCH. Null model: AR(',ar.ord,')- GARCH(',arch.ord,',',garch.ord,')',sep='')
  }else{
    METHOD <- paste('supLM test ARMA vs TARMA with GARCH innovations. Null model: ARMA(',ar.ord,',',ma.ord,')- GARCH(',arch.ord,',',garch.ord,')',sep='')
  }
  structure(list(statistic=test.stat, parameter=c(thd), test.v=test.v+shift,
                 thd.range=thd.range+shift, fit.ARMA=fit, fit.GARCH=fit.g, sigma2=s2,
                  data.name=DNAME,prop=mean(Il), p.value=NULL,
                 method=METHOD,d=d,pa=pa, dfree=dfree),class=c('TARMAtest','htest'))
}
## ****************************************************************************
