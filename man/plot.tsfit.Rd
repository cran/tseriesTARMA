% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.tsfit.R
\name{plot.tsfit}
\alias{plot.tsfit}
\title{Plot from fitted/forecasted time series models.}
\usage{
\method{plot}{tsfit}(
  x,
  fit,
  plot.fit = TRUE,
  fore = NULL,
  lcols = list(series = "lightblue", fit = 4, pred = "red4", band = "red"),
  ptype = list(series = 20, fit = 1, pred = 16, band = 20),
  ltype = list(series = 1, fit = 1, pred = 1, band = 1),
  lwdt = 2,
  ...
)
}
\arguments{
\item{x}{A time series}

\item{fit}{A time series model fitted upon \code{x}, e.g. an object obtained with \code{TARMA.fit}}

\item{plot.fit}{Logical. If \code{TRUE} adds the fitted values from the model.}

\item{fore}{Forecast derived from \code{fit}, e.g. an object obtained with \code{predict.TARMA}
if not null it adds the prediction together with its confidence bands.}

\item{lcols}{List of colours for the plots.}

\item{ptype}{List of point types (\code{pch}) for the plots.}

\item{ltype}{List of line types (\code{lty}) for the plots.}

\item{lwdt}{A common line width for the plots.}

\item{\dots}{Additional graphical parameters.}
}
\value{
No return value, called for side effects
}
\description{
Plots a time series model fit, possibly together with its forecast and forecast bands.
}
\examples{
## a TARMA(1,1,1,1) model
set.seed(13)
x    <- TARMA.sim(n=200, phi1=c(0.5,-0.5), phi2=c(0.0,0.5), theta1=-0.5, theta2=0.7, d=1, thd=0.2)
fit1 <- TARMA.fit(x,tar1.lags = 1, tar2.lags = 1, tma1.lags = 1, tma2.lags = 1, d=1,
         estimate.thd = FALSE, threshold=0.2)
xp1  <- predict(fit1,x,n.ahead=5)

# plots both the fitted and the forecast
plot.tsfit(x,fit=fit1,fore=xp1);grid();

# plots only the forecast
plot.tsfit(x,fit=fit1,plot.fit=FALSE,fore=xp1);grid();

# plots only the fitted
plot.tsfit(x,fit=fit1);grid();

}
\seealso{
\code{\link{TARMA.fit}} and \code{\link{TARMA.fit2}} for TARMA modelling.
\code{\link{predict.TARMA}} for prediction and forecasting.
}
\author{
Simone Giannerini, \email{simone.giannerini@unibo.it}

Greta Goracci, \email{greta.goracci@unibz.it}
}
