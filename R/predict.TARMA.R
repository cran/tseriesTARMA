#' Forecast from fitted TARMA models.
#'
#'  Forecasting with TARMA models
#'
#' @param object   A \code{TARMA} fit upon x.
#' @param x     The fitted time series.
#' @param n.ahead The number of steps ahead for which prediction is required.
#' @param n.sim   The number of Monte Carlo replications used to simulate the prediction density.
#' @param pred.matrix Logical. if \code{TRUE} prints also the whole simulated prediction density for each prediction horizon from \code{1} to \code{n.ahead}.
#' @param quant   Vector of quantiles (in the interval \code{[0, 1]}) to be computed upon the prediction density.
#' @param \dots   Additional arguments.
#'
#' @details If \code{n.ahead = 0} it gives the fitted values from the model.
#' If the fit is from \code{TARMA.fit2} and includes covariates, these are ignored.
#' @return
#' A list with components `pred.matrix`, `pred`, and `pred.interval`. The latter two are `ts` objects that contain the prediction and the quantiles of the prediction density, respectively.
#' If `pred.matrix = TRUE` then the prediction density from which the quantiles are computed is also returned. 
#' @author Simone Giannerini, \email{simone.giannerini@@unibo.it}
#' @author Greta Goracci, \email{greta.goracci@@unibo.it}
#' @references
#' * \insertRef{Gia21}{tseriesTARMA}
#'
#' @importFrom stats deltat
#' @method predict TARMA
#' @export
#' @seealso \code{\link{TARMA.fit}} and \code{\link{TARMA.fit2}} for TARMA modelling. \code{\link{plot.tsfit}} for plotting TARMA fits and forecasts.
#' @examples
#' ## a TARMA(1,1,1,1) model
#' set.seed(13)
#' x1   <- TARMA.sim(n=200, phi1=c(0.5,-0.5), phi2=c(0.0,0.5), theta1=-0.5, theta2=0.7, d=1, thd=0.2)
#' fit1 <- TARMA.fit(x1, method='L-BFGS-B',tar1.lags = 1, tar2.lags = 1, tma1.lags = 1, 
#'         tma2.lags = 1, d=1, threshold=0.2)
#' xp1  <- predict(fit1,x1,n.ahead=2)
#' xp1
## ***************************************************************************

predict.TARMA <- function(object,x,n.ahead=0, n.sim=1e3, quant=c(0.05,0.95), pred.matrix=FALSE,...){
    n    <- length(x)
    xtsp <- tsp(x)
    neff <- object$fit$nobs # effective sample size
    d    <- object$d
    k    <- max(object$arlag,object$tlag1,object$tlag2,object$mlag1,object$mlag2,d) # number of discarded observations
    eps  <- object$fit$residuals
    sdx  <- sqrt(object$fit$sigma2)
    if(n.ahead==0){
        fit1  <- x[-(1:k)]- eps
        res <- NULL
        xp  <- NULL
    }else{
        nq     <- length(quant)
        ma.ord <- max(object$mlag1,object$mlag2)
        thd    <- object$thd
        if(x[n+1-d]<=thd){
            lag.g <- object$tlag1
        } else {lag.g <- object$tlag2}
        phi1   <- object$phi1
        phi2   <- object$phi2
        theta1 <- object$theta1
        theta2 <- object$theta2
        res <- replicate(n=n.sim,expr=TARMA.sim(n=n.ahead,rand.gen=function(n,...) {rnorm(n,sd=sdx)},
         n.start=0,xstart=list(ar=(x[n-lag.g+1]), ma=eps[neff:(neff-ma.ord+1)])
        ,phi1=phi1, phi2=phi2, theta1=theta1, theta2=theta2, d=d, thd=thd, s1=1,s2=1))
        res  <- matrix(res,n.ahead,n.sim)
        xp   <- ts(t(apply(res,MARGIN=1,FUN=quantile,probs=quant)),start=xtsp[2L]+deltat(x), frequency=xtsp[3L])
        fit1 <- ts(c(t(apply(res,MARGIN=1,FUN=mean))),start=xtsp[2L]+deltat(x), frequency=xtsp[3L])
        if(!pred.matrix){
            res <- NULL
        }
    }
    return(list(pred.matrix=res, pred = fit1, pred.interval = xp))
}
