% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TARMA.fit.R
\name{TARMA.fit}
\alias{TARMA.fit}
\title{TARMA Modelling of Time Series}
\usage{
TARMA.fit(
  x,
  tar1.lags = c(1),
  tar2.lags = c(1),
  tma1.lags = c(1),
  tma2.lags = c(1),
  estimate.thd = TRUE,
  threshold,
  d = 1,
  pa = 0.25,
  pb = 0.75,
  method = c("L-BFGS-B", "solnp", "lbfgsb3c", "robust"),
  alpha = 0,
  qu = c(0.05, 0.95),
  optim.control = list(trace = 0),
  ...
)
}
\arguments{
\item{x}{A univariate time series.}

\item{tar1.lags}{Vector of AR lags for the lower regime. It can be a subset of \code{1 ... p1 = max(tar1.lags)}.}

\item{tar2.lags}{Vector of AR lags for the upper regime. It can be a subset of \code{1 ... p2 = max(tar2.lags)}.}

\item{tma1.lags}{Vector of MA lags for the lower regime. It can be a subset of \code{1 ... q1 = max(tma1.lags)}.}

\item{tma2.lags}{Vector of MA lags for the upper regime. It can be a subset of \code{1 ... q2 = max(tma2.lags)}.}

\item{estimate.thd}{Logical. If \code{TRUE} estimates the threshold over the threshold range specified by \code{pa} and \code{pb}.}

\item{threshold}{Threshold parameter. Used only if \code{estimate.thd = FALSE}.}

\item{d}{Delay parameter. Defaults to \code{1}.}

\item{pa}{Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
The default is \code{0.25}}

\item{pb}{Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
The default is \code{0.75}}

\item{method}{Optimization/fitting method, can be one of \cr
\code{"L-BFGS-B", "solnp", "lbfgsb3c", "robust"}.}

\item{alpha}{Real positive number. Tuning parameter for robust estimation. Only used if \code{method} is \code{"robust"}.}

\item{qu}{Quantiles for initial trimmed estimation. Tuning parameter for robust estimation. Only used if \code{method} is \code{"robust"}.}

\item{optim.control}{List of control parameters for the optimization method.}

\item{\dots}{Additional arguments for the optimization.}
}
\value{
A list of class \code{TARMA} with components:
\itemize{
\item \code{fit} - List with the following components \cr
\itemize{
\item \code{coef} - Vector of estimated parameters which can be extracted by the coef method.
\item \code{sigma2} - Estimated innovation variance.
\item \code{var.coef} - The estimated variance matrix of the coefficients coef, which can be extracted by the vcov method
\item \code{residuals} - Vector of residuals from the fit.
\item \code{nobs} - Effective sample size used for fitting the model.
}
\item \code{se} - Standard errors for the parameters. Note that they are computed conditionally \cr
upon the threshold so that they are generally smaller than the true ones.
\item \code{thd}    - Estimated threshold.
\item \code{aic}    - Value of the AIC for the minimised least squares criterion over the threshold range.
\item \code{bic}    - Value of the BIC for the minimised least squares criterion over the threshold range.
\item \code{rss}    - Minimised value of the target function. Coincides with the residual sum of squares for ordinary least squares estimation.
\item \code{rss.v}  - Vector of values of the rss over the threshold range.
\item \code{thd.range} - Vector of values of the threshold range.
\item \code{d}      - Delay parameter.
\item \code{phi1}   - Estimated AR parameters for the lower regime.
\item \code{phi2}   - Estimated AR parameters for the upper regime.
\item \code{theta1} - Estimated MA parameters for the lower regime.
\item \code{theta2} - Estimated MA parameters for the upper regime.
\item \code{tlag1}  - TAR lags for the lower regime
\item \code{tlag2}  - TAR lags for the upper regime
\item \code{mlag1}  - TMA lags for the lower regime
\item \code{mlag2}  - TMA lags for the upper regime
\item \code{method} - Estimation method.
\item \code{alpha}  - Tuning parameter for robust estimation.
\item \code{qu}     - Tuning parameter for robust estimation.
\item \code{call}   - The matched call.
\item \code{convergence} - Convergence code from the optimization routine.
}
}
\description{
\loadmathjax
Implements a Least Squares fit of full subset two-regime \code{TARMA(p1,p2,q1,q2)} model to a univariate time series
}
\details{
Implements the Least Squares fit of the following two-regime \code{TARMA(p1,p2,q1,q2)} process: \cr
\mjdeqn{X_{t} = \left\lbrace
    \begin{array}{ll}
\phi_{1,0} + \sum_{i \in I_1} \phi_{1,i} X_{t-i} + \sum_{j \in M_1} \theta_{1,j} \varepsilon_{t-j} + \varepsilon_{t} & \mathrm{if } X_{t-d} \leq \mathrm{thd} \\\\\\
 &\\\\\\
\phi_{2,0} + \sum_{i \in I_2} \phi_{2,i} X_{t-i} + \sum_{j \in M_2} \theta_{2,j} \varepsilon_{t-j} + \varepsilon_{t} & \mathrm{if } X_{t-d} > \mathrm{thd}
\end{array}
\right. }{X[t] =
 \phi[1,0] + \Sigma_{i in I_1} \phi[1,i] X[t-i] + \Sigma_{j in M_1} \theta[1,j] \epsilon[t-j] + \epsilon[t] --  if X[t-d] <= thd
 \phi[2,0] + \Sigma_{i in I_2} \phi[2,i] X[t-i] + \Sigma_{j in M_2} \theta[2,j] \epsilon[t-j] + \epsilon[t] --  if X[t-d] > thd}
where  \mjeqn{\phi_{1,i}}{\phi[1,i]} and \mjeqn{\phi_{2,i}}{\phi[2,i]} are the TAR parameters for the lower and upper regime, respectively, and
\code{I1 = tar1.lags} and \code{I2 = tar2.lags} are the corresponding vectors of TAR lags.
\mjeqn{\theta_{1,j}}{\theta[1,j]} and \mjeqn{\theta_{2,j}}{\theta[2,j]} are the TMA parameters
and \mjeqn{j \in M_1, M_2}{j in M_1, M_2}, where \code{M1 = tma1.lags} and \code{M2 = tma2.lags}, are the vectors of TMA lags. \cr
The most demanding routines have been reimplemented in Fortran and dynamically loaded.
}
\section{Fitting methods}{

\code{method} has the following options: \cr
\describe{
\item{\code{L-BFGS-B}}{Calls the corresponding method of \code{\link{optim}}. Linear ergodicity constraints are imposed.}
\item{\code{solnp}}{Calls the function \code{\link[Rsolnp]{solnp}}. It is a nonlinear optimization using augmented Lagrange method with
linear and nonlinear inequality bounds. This allows to impose all the ergodicity constraints so that in theory it always
return an ergodic solution. In practice the solution should be checked since this is a local solver and there is no guarantee
that the minimum has been reached.}
\item{\code{lbfgsb3c}}{Calls the function \code{\link[lbfgsb3c]{lbfgsb3c}} in package \code{lbfgsb3c}. Improved version of the L-BFGS-B in \code{\link{optim}}.}
\item{\code{robust}}{Robust M-estimator of Ferrari and La Vecchia \insertCite{Fer12}{tseriesTARMA}. Based on the L-BFGS-B in \code{\link{optim}} and an additional iterative reweighted least squares step to estimate the robust weights.
Uses the tuning parameters \code{alpha} and \code{qu}. Robust standard errors are derived from the sandwich estimator of the variance/covariance matrix of the estimates}
}
Where possible, the ergodicity bounds are imposed to the optimization routines but there is no guarantee that the solution will be ergodic so that it is
advisable to check the fitted parameters.
}

\examples{
\donttest{
## a TARMA(1,1,1,1) model
set.seed(13)
x    <- TARMA.sim(n=200, phi1=c(0.5,-0.5), phi2=c(0.0,0.5), theta1=-0.5, theta2=0.7, d=1, thd=0.2)
fit1 <- TARMA.fit(x,tar1.lags=1, tar2.lags=1, tma1.lags=1, tma2.lags=1, d=1)
}
## --------------------------------------------------------------------------
## In the following examples the threshold is fixed to speed up computations
## --------------------------------------------------------------------------

## --------------------------------------------------------------------------
## Least Squares fit
## --------------------------------------------------------------------------

set.seed(26)
n    <- 200
y    <- TARMA.sim(n=n, phi1=c(0.6,0.6), phi2=c(-1.0,0.4), theta1=-0.7, theta2=0.5, d=1, thd=0.2)

fit1 <- TARMA.fit(y,tar1.lags=1, tar2.lags=1, tma1.lags=1, tma2.lags=1, d=1, 
       estimate.thd=FALSE, threshold=0.2)
fit1

## ---------------------------------------------------------------------------
## Contaminating the data with one additive outlier
## ---------------------------------------------------------------------------
x     <- y           # contaminated series
x[54] <- x[54] + 10

## ---------------------------------------------------------------------------
## Compare the non-robust LS fit with the robust fit
## ---------------------------------------------------------------------------

fitls  <- TARMA.fit(x,tar1.lags=1, tar2.lags=1, tma1.lags=1, tma2.lags=1, d=1,
            estimate.thd=FALSE, threshold=0.2)
fitrob <- TARMA.fit(x,tar1.lags=1, tar2.lags=1, tma1.lags=1, tma2.lags=1, d=1,
            method='robust',alpha=0.7,qu=c(0.1,0.95),estimate.thd = FALSE, threshold=0.2)

par.true <- c(0.6,0.6,-1,0.4,-0.7,0.5)
pnames   <- c("int.1", "ar1.1", "int.2", "ar2.1", "ma1.1", "ma2.1")
names(par.true) <- pnames

par.ls  <- round(fitls$fit$coef,2)  # Least Squares
par.rob <- round(fitrob$fit$coef,2) # robust

rbind(par.true,par.ls,par.rob)
}
\references{
\itemize{
\item \insertRef{Gia21}{tseriesTARMA}
\item \insertRef{Cha19}{tseriesTARMA}
\item \insertRef{Gor23b}{tseriesTARMA}
\item \insertRef{Fer12}{tseriesTARMA}
}
}
\seealso{
\code{\link{TARMA.fit2}} for Maximum Likelihood estimation of TARMA
models with common MA part. \code{\link{print.TARMA}} for print methods for \code{TARMA} fits.
\code{\link{predict.TARMA}} for prediction and forecasting. \code{\link{plot.tsfit}} for plotting TARMA fits and forecasts.
}
\author{
Simone Giannerini, \email{simone.giannerini@unibo.it}

Greta Goracci, \email{greta.goracci@unibz.it}
}
