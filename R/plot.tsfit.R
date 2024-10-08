#' Plot from fitted/forecasted time series models.
#'
#'  Plots a time series model fit, possibly together with its forecast and forecast bands.
#'
#' @param x        A time series
#' @param fit      A time series model fitted upon \code{x}, e.g. an object obtained with \code{TARMA.fit}
#' @param plot.fit Logical. If \code{TRUE} adds the fitted values from the model.
#' @param fore     Forecast derived from \code{fit}, e.g. an object obtained with \code{predict.TARMA}
#' if not null it adds the prediction together with its confidence bands.
#' @param lcols    List of colours for the plots.
#' @param ptype    List of point types (\code{pch}) for the plots.
#' @param ltype    List of line types (\code{lty}) for the plots.
#' @param lwdt     A common line width for the plots.
#' @param \dots    Additional graphical parameters.
#' @return
#'   No return value, called for side effects
#' @author Simone Giannerini, \email{simone.giannerini@@uniud.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @importFrom stats ts.plot predict
#' @rawNamespace export(plot.tsfit)
#' @export
#' @seealso \code{\link{TARMA.fit}} and \code{\link{TARMA.fit2}} for TARMA modelling.
#' \code{\link{predict.TARMA}} for prediction and forecasting.
#' @examples
#' ## a TARMA(1,1,1,1) model
#' set.seed(13)
#' x    <- TARMA.sim(n=200, phi1=c(0.5,-0.5), phi2=c(0.0,0.5), theta1=-0.5, theta2=0.7, d=1, thd=0.2)
#' fit1 <- TARMA.fit(x,tar1.lags = 1, tar2.lags = 1, tma1.lags = 1, tma2.lags = 1, d=1, threshold=0.2)
#' xp1  <- predict(fit1,x,n.ahead=5)
#'
#' # plots both the fitted and the forecast
#' plot.tsfit(x,fit=fit1,fore=xp1);grid();
#'
#' # plots only the forecast
#' plot.tsfit(x,fit=fit1,plot.fit=FALSE,fore=xp1);grid();
#'
#' # plots only the fitted
#' plot.tsfit(x,fit=fit1);grid();
#' 
## '***************************************************************************

plot.tsfit <- function(x,fit, plot.fit=TRUE, fore=NULL, lcols=list(series='lightblue',fit=4,pred='red4',band='red'),
ptype=list(series=20,fit=1,pred=16,band=20), ltype=list(series=1,fit=1,pred=1,band=1), lwdt=2, ...){
    if(inherits(fit, "TARMA")){
        if(plot.fit){
            xh <- predict(fit,x)$pred
            y  <- ts.intersect(x,xh)
            lcol <- c(lcols$series,lcols$fit)
            ptyp <- c(ptype$series,ptype$fit)
            ltyp <- c(ltype$series,ltype$fit)
        }else{
            y    <- x
            lcol <- c(lcols$series)
            ptyp <- ptype$series
            ltyp <- ltype$series
        }
        ny  <- NCOL(y)
        if(is.null(fore)){
            ts.plot(y,gpars=list(col=lcol,lwd=lwdt,pch=ptyp,lty=ltyp),...);
        }else{
            npred <- length(fore$pred)
            lcol <- c(lcol,lcols$pred,lcols$band,lcols$band)
            ptyp <- c(ptyp,ptype$pred,ptype$band,ptype$band)
            ltyp <- c(ltyp,ltype$pred,ltype$band,ltype$band)
            ts.plot(y,fore$pred,fore$pred.interval,gpars=list(col=lcol,lwd=lwdt,pch=ptyp,lty=ltyp),...);
        }
    }else if(inherits(fit,'Arima')){
       stop('method not yet implemented for Arima objects')
    }
return(invisible())
}
