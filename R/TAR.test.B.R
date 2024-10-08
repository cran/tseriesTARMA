#'  AR versus TAR bootstrap supLM test for nonlinearity
#'
#'  Implements various bootstrap supremum Lagrange Multiplier tests for a AR specification versus
#'  a TAR specification.
#'
#' @param x         A univariate time series.
#' @param B         Integer. Number of bootstrap resamples. Defaults to 1000.
#' @param pa        Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
#'                  The default is \code{0.25}
#' @param pb        Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
#'                  The default is \code{0.75}
#' @param ar.ord    Order of the AR part.
#' @param d         Delay parameter. Defaults to \code{1}.
#' @param btype     Bootstrap type, can be one of \code{'iid','wb.h','wb.r','wb.n'}, see Details.
#' @param \dots     Additional arguments to be passed to \code{arima}.
#'
#' @details
#' Implements the bootstrap version of \code{\link{TAR.test}} the supremum Lagrange Multiplier test to test an AR specification versus a TARMA specification.
#' The option \code{btype} specifies the type of bootstrap as follows:
#'  \describe{
#'  \item{\code{iid}}{Residual iid bootstrap. See \insertCite{Gia22}{tseriesTARMA}, \insertCite{Gia23}{tseriesTARMA}.}
#'  \item{\code{wb.h}}{Stochastic permutation of \insertCite{Han96}{tseriesTARMA}.}
#'  \item{\code{wb.r}}{Residual wild bootstrap with Rademacher auxiliary distribution. See \insertCite{Gia22}{tseriesTARMA}, \insertCite{Gia23}{tseriesTARMA}.}
#'  \item{\code{wb.n}}{Residual wild bootstrap with Normal auxiliary distribution. See \insertCite{Gia22}{tseriesTARMA}, \insertCite{Gia23}{tseriesTARMA}.}
#' }
#'
#' @return
#'   A list of class \code{htest} with components:
#' \describe{
#'  \item{\code{statistic}}{The value of the supLM statistic.}
#'  \item{\code{parameter}}{A named vector: \code{threshold} is the value that maximises the Lagrange Multiplier values.}
#'  \item{\code{test.v}}{Vector of values of the LM statistic for each threshold given in \code{thd.range}.}
#'  \item{\code{thd.range}}{Range of values of the threshold.}
#'  \item{\code{fit}}{The null model: AR fit over \code{x}.}
#'  \item{\code{sigma2}}{Estimated innovation variance from the AR fit.}
#'  \item{\code{data.name}}{A character string giving the name of the data.}
#'  \item{\code{prop}}{Proportion of values of the series that fall in the lower regime.}
#'  \item{\code{p.value}}{The bootstrap p-value of the test.}
#'  \item{\code{method}}{A character string indicating the type of test performed.}
#'  \item{\code{Tb}}{The bootstrap null distribution.}
#'}
#' @importFrom stats coef residuals
#' @export
#' @author Simone Giannerini, \email{simone.giannerini@@uniud.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Gia22}{tseriesTARMA}
#' * \insertRef{Gia23}{tseriesTARMA}
#' * \insertRef{Gor23}{tseriesTARMA}
#' * \insertRef{Gia21}{tseriesTARMA}
#' * \insertRef{Han96}{tseriesTARMA}
#'
#' @seealso \code{\link{TAR.test}} for the heteroskedastic robust asymptotic test.  \code{\link{TARMAGARCH.test}} for the
#' robust version of the test with respect to GARCH innovations. \code{\link{TARMA.sim}} to simulate from a TARMA process.
#'
#' @examples
#' ## a TAR(1,1) where the threshold effect is on the AR parameters
#' set.seed(123)
#' x1 <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0.0,0.8), theta1=0, theta2=0, d=1, thd=0.2)
#' TAR.test.B(x1, ar.ord=1, d=1)
#' TAR.test.B(x1, ar.ord=1, d=1, btype='wb.r')
#' TAR.test.B(x1, ar.ord=1, d=1, btype='wb.h')
#' 
#' ## a AR(1)
#' x2 <- arima.sim(n=100, model=list(order = c(1,0,0),ar=0.5))
#' TAR.test.B(x2, ar.ord=1, d=1)
#' TAR.test.B(x2, ar.ord=1, d=1, btype='wb.r')
#' TAR.test.B(x2, ar.ord=1, d=1, btype='wb.h')
#' 
## '***************************************************************************

TAR.test.B <- function(x, B=1000 ,pa=.25, pb=.75, ar.ord, d=1, 
                       btype = c('iid','wb.h','wb.r','wb.n'),...){
  
  btype  <- match.arg(btype)
  DNAME  <- deparse(substitute(x))
  n      <- length(x)
  test.a <- TAR.test(x,pa=pa, pb=pb, ar.ord, d)
  test.stat <- test.a$statistic
  resi      <- scale(test.a$fit$residuals,center=TRUE,scale=FALSE) # centered residuals for the bootstrap
  int.p     <- test.a$fit$coefficients['Intercept']
  ar.p      <- test.a$fit$coefficients[-1]
  
  k    <- max(ar.ord,d) # number of discarded obs
  neff <- n - k
  a    <- ceiling((neff-1)*pa)
  b    <- floor((neff-1)*pb)
  nr   <- b-a+1
  
  # Bootstrap
  # **************************************************************************
  n.start <- floor(n/3)
  np     <- neff + n.start # transient
  if(btype == 'wb.r'){
    eta <- matrix(sample(c(-1,1),replace=TRUE, size=np*B),np,B)
    resi.b    <- c(resi[1:n.start],resi)*eta
    smeth     <- 'wild bootstrap (Rademacher)'
  }else
    if(btype == 'wb.n'){
      eta <- matrix(rnorm(np*B),np,B)
      resi.b    <- c(resi[1:n.start],resi)*eta
      smeth     <- 'wild bootstrap (Gaussian)'
    }else
      if(btype == 'iid'){
        resi.b    <- matrix(sample(resi,size=np*B,replace=TRUE),np,B)
        smeth     <- 'i.i.d. bootstrap'
      }
  # ******************************************************************************
  if(btype == 'wb.h'){
    smeth     <- 'stochastic perturbation (Hansen)'
    #       SUBROUTINE ARvsTAR_HB(x,n,d,p,a,b,neff,nrep,test,testb)
    dum <- matrix(0,B,2)
    storage.mode(dum) <- 'double'
    Tb   <- .Fortran('arvstar_hb',as.double(x),as.integer(n),as.integer(d),
                     as.integer(ar.ord),as.integer(a),as.integer(b),as.integer(neff),
                     as.integer(B),test=double(2),testb = dum,PACKAGE='tseriesTARMA')$testb
  }else{
    x.b <- resi.b
    for(i in (ar.ord+1):np){
      x.b[i,] <- int.p + c(t(x.b[(i-1):(i-ar.ord),,drop=FALSE])%*%ar.p) + resi.b[i,] #
    }
    x.b  <- ts(x.b[-(1:(n.start-k)),]) # bootstrap resamples [n x B]
    Tb   <- matrix(NA,B,2)
    for(k in 1:B){
      xb      <- x.b[,k]
      res     <- TAR.test(xb,pa=pa, pb=pb, ar.ord, d)
      Tb[k, ] <- res$statistic
    }
  }
  dfree        <- 1 + ar.ord
  colnames(Tb) <- c('sLM','sLMh')
  p.value      <- rowMeans(t(Tb)>test.stat,na.rm=TRUE)
  METHOD       <- paste(smeth,' supLM test AR vs TAR. Null model: AR(',ar.ord,')',sep='')
  structure(list(test.v=test.a$test.v,
                 thd.range=test.a$thd.range, fit=test.a$fit, sigma2=test.a$sigma2,
                 parameter=test.a$parameter, data.name=DNAME, prop=test.a$prop, statistic=test.a$statistic,
                 p.value=p.value[1], method=METHOD,Tb=Tb,pval=p.value,
                 d=d,pa=pa, dfree=dfree),class=c('TARMAtest','htest'))
}
## ****************************************************************************









