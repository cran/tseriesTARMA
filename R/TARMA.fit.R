#'  TARMA Modelling of Time Series
#'
#' @description \loadmathjax
#'  Implements a Least Squares fit of full subset two-regime \code{TARMA(p1,p2,q1,q2)} model to a univariate time series
#'
#' @param x   A univariate time series.
#' @param tar1.lags Vector of AR lags for the lower regime. It can be a subset of \code{1 ... p1 = max(tar1.lags)}.
#' @param tar2.lags Vector of AR lags for the upper regime. It can be a subset of \code{1 ... p2 = max(tar2.lags)}.
#' @param tma1.lags Vector of MA lags for the lower regime. It can be a subset of \code{1 ... q1 = max(tma1.lags)}.
#' @param tma2.lags Vector of MA lags for the upper regime. It can be a subset of \code{1 ... q2 = max(tma2.lags)}.
#' @param threshold Threshold parameter. If \code{NULL} estimates the threshold over the threshold range specified by 
#' \code{pa} and \code{pb}.
#' @param d  Delay parameter. Defaults to \code{1}.
#' @param pa  Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
#' The default is \code{0.25}
#' @param pb  Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
#' The default is \code{0.75}
#' @param method  Optimization/fitting method, can be one of \code{"L-BFGS-B", "solnp", "lbfgsb3c", "robust", "trimmed"}.
#' @param alpha  Real positive number. Tuning parameter for robust estimation. Only used if \code{method} is \code{robust}.
#' @param qu     Quantiles for (initial) trimmed estimation. Tuning parameter for robust estimation. Only used if \code{method} is either \code{robust} or \code{trimmed}.
#' @param innov  Innovation density for robust estimation. can be one of \code{"norm", "student"}. Only used if \code{method} is \code{"robust"}.
#' @param optim.control List of control parameters for the main optimization method.
#' @param irls.control List of control parameters for the irls optimization method (see details).
#' @param \dots  Additional arguments.
#'
#' @details
#'  Implements the Least Squares fit of the following two-regime \code{TARMA(p1,p2,q1,q2)} process: \cr
#' \mjdeqn{X_{t} = \left\lbrace
#'     \begin{array}{ll}
#' \phi_{1,0} + \sum_{i \in I_1} \phi_{1,i} X_{t-i} + \sum_{j \in M_1} \theta_{1,j} \varepsilon_{t-j} + \varepsilon_{t} & \mathrm{if } X_{t-d} \leq \mathrm{thd} \\\\\\
#'  &\\\\\\
#' \phi_{2,0} + \sum_{i \in I_2} \phi_{2,i} X_{t-i} + \sum_{j \in M_2} \theta_{2,j} \varepsilon_{t-j} + \varepsilon_{t} & \mathrm{if } X_{t-d} > \mathrm{thd}
#' \end{array}
#' \right. }{X[t] =
#'  \phi[1,0] + \Sigma_{i in I_1} \phi[1,i] X[t-i] + \Sigma_{j in M_1} \theta[1,j] \epsilon[t-j] + \epsilon[t] --  if X[t-d] <= thd
#'  \phi[2,0] + \Sigma_{i in I_2} \phi[2,i] X[t-i] + \Sigma_{j in M_2} \theta[2,j] \epsilon[t-j] + \epsilon[t] --  if X[t-d] > thd}
#' where  \mjeqn{\phi_{1,i}}{\phi[1,i]} and \mjeqn{\phi_{2,i}}{\phi[2,i]} are the TAR parameters for the lower and upper regime, respectively, and
#' \code{I1 = tar1.lags} and \code{I2 = tar2.lags} are the corresponding vectors of TAR lags.
#' \mjeqn{\theta_{1,j}}{\theta[1,j]} and \mjeqn{\theta_{2,j}}{\theta[2,j]} are the TMA parameters
#' and \mjeqn{j \in M_1, M_2}{j in M_1, M_2}, where \code{M1 = tma1.lags} and \code{M2 = tma2.lags}, are the vectors of TMA lags. \cr
#' The most demanding routines have been reimplemented in Fortran and dynamically loaded.
#' @section Fitting methods:
#' \code{method} has the following options: \cr
#' \describe{
#'  \item{\code{L-BFGS-B}}{Calls the corresponding method of \code{\link{optim}}. Linear ergodicity constraints are imposed.}
#'  \item{\code{solnp}}{Calls the function \code{\link[Rsolnp]{solnp}}. It is a nonlinear optimization using augmented Lagrange method with
#'     linear and nonlinear inequality bounds. This allows to impose all the ergodicity constraints so that in theory it always
#'     return an ergodic solution. In practice the solution should be checked since this is a local solver and there is no guarantee
#'     that the minimum has been reached.}
#'  \item{\code{lbfgsb3c}}{Calls the function \code{\link[lbfgsb3c]{lbfgsb3c}} in package \code{lbfgsb3c}. Improved version of the L-BFGS-B in \code{\link{optim}}.}
#'  \item{\code{robust}}{Robust M-estimator of Ferrari and La Vecchia \insertCite{Fer12}{tseriesTARMA}. Based on the L-BFGS-B in \code{\link{optim}} and an additional iterative re-weighted least squares step to estimate the robust weights.
#'        Uses the tuning parameters \code{alpha} and \code{qu}. Robust standard errors are derived from the sandwich estimator of the variance/covariance matrix of the estimates.
#'        The IRLS step can be controlled through the parameters \code{maxiter} (maximum number of iterations) and \code{tol} (target tolerance). These can be passed using \code{irls.control}.}
#'  \item{\code{trimmed}}{Experimental: Estimator based on trimming the sample using the tuning parameters \code{qu} (lower and upper quantile).}
#' }
#' Where possible, the conditions for ergodicity and invertibility are imposed to the optimization routines but there is no guarantee that the solution will be ergodic and invertible so that it is
#' advisable to check the fitted parameters.
#' @return
#'   A list of class \code{TARMA} with components:
#' \itemize{
#' \item \code{fit} - List with the following components \cr
#'     \itemize{
#'       \item \code{coef} - Vector of estimated parameters which can be extracted by the coef method.
#'       \item \code{sigma2} - Estimated innovation variance.
#'       \item \code{var.coef} - The estimated variance matrix of the coefficients coef, which can be extracted by the vcov method
#'       \item \code{residuals} - Vector of residuals from the fit.
#'       \item \code{nobs} - Effective sample size used for fitting the model.
#'     }
#'  \item \code{se} - Standard errors for the parameters. Note that they are computed conditionally upon the threshold so that they are generally smaller than the true ones.
#'  \item \code{thd}    - Estimated threshold.
#'  \item \code{aic}    - Value of the AIC for the minimised least squares criterion over the threshold range.
#'  \item \code{bic}    - Value of the BIC for the minimised least squares criterion over the threshold range.
#'  \item \code{rss}    - Minimised value of the target function. Coincides with the residual sum of squares for ordinary least squares estimation. 
#'  \item \code{rss.v}  - Vector of values of the rss over the threshold range.
#'  \item \code{thd.range} - Vector of values of the threshold range.
#'  \item \code{d}      - Delay parameter.
#'  \item \code{phi1}   - Estimated AR parameters for the lower regime.
#'  \item \code{phi2}   - Estimated AR parameters for the upper regime.
#'  \item \code{theta1} - Estimated MA parameters for the lower regime.
#'  \item \code{theta2} - Estimated MA parameters for the upper regime.
#'  \item \code{tlag1}  - TAR lags for the lower regime
#'  \item \code{tlag2}  - TAR lags for the upper regime
#'  \item \code{mlag1}  - TMA lags for the lower regime
#'  \item \code{mlag2}  - TMA lags for the upper regime
#'  \item \code{method} - Estimation method.
#'  \item \code{innov}  - Innovation density model.
#'  \item \code{alpha}  - Tuning parameter for robust estimation.
#'  \item \code{qu}     - Tuning parameter for robust estimation.
#'  \item \code{call}   - The matched call.
#'  \item \code{convergence} - Convergence code from the optimization routine.
#'  \item \code{innovpar} - Parameter vector for the innovation density. Defaults to \code{NULL}.
#'}
#' @author Simone Giannerini, \email{simone.giannerini@@uniud.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Gia21}{tseriesTARMA}
#' * \insertRef{Cha19}{tseriesTARMA}
#' * \insertRef{Gor23b}{tseriesTARMA}
#' * \insertRef{Fer12}{tseriesTARMA}
#'
#' @importFrom zoo is.zoo zoo index
#' @importFrom methods new
#' @importFrom Rsolnp solnp
#' @importFrom lbfgsb3c lbfgsb3c
#' @importFrom stats is.ts window median quantile end frequency tsp tsp<- pnorm var
#' @importFrom Matrix forceSymmetric
#' @importFrom fitdistrplus fitdist
#' @export
#' @seealso \code{\link{TARMA.fit2}} for Maximum Likelihood estimation of TARMA
#' models with common MA part. \code{\link{print.TARMA}} for print methods for \code{TARMA} fits.
#' \code{\link{predict.TARMA}} for prediction and forecasting. \code{\link{plot.tsfit}} for plotting TARMA fits and forecasts.
#' @examples
#' \donttest{
#' ## a TARMA(1,1,1,1) model
#' set.seed(13)
#' x    <- TARMA.sim(n=200, phi1=c(0.5,-0.5), phi2=c(0.0,0.5), theta1=-0.5, theta2=0.7, d=1, thd=0.2)
#' fit1 <- TARMA.fit(x,tar1.lags=1, tar2.lags=1, tma1.lags=1, tma2.lags=1, d=1)
#' }
#' ## --------------------------------------------------------------------------
#' ## In the following examples the threshold is fixed to speed up computations
#' ## --------------------------------------------------------------------------
#' 
#' ## --------------------------------------------------------------------------
#' ## Least Squares fit
#' ## --------------------------------------------------------------------------
#'
#' set.seed(26)
#' n    <- 200
#' y    <- TARMA.sim(n=n, phi1=c(0.6,0.6), phi2=c(-1.0,0.4), theta1=-0.7, theta2=0.5, d=1, thd=0.2)
#'
#' fit1 <- TARMA.fit(y,tar1.lags=1, tar2.lags=1, tma1.lags=1, tma2.lags=1, d=1, threshold=0.2)
#' fit1
#' 
#' ## ---------------------------------------------------------------------------
#' ## Contaminate the data with one additive outlier
#' ## ---------------------------------------------------------------------------
#' x     <- y           # contaminated series
#' x[54] <- x[54] + 10
#' 
#' ## ---------------------------------------------------------------------------
#' ## Compare the non-robust LS fit with the robust fit
#' ## ---------------------------------------------------------------------------
#'
#' fitls  <- TARMA.fit(x,tar1.lags=1, tar2.lags=1, tma1.lags=1, tma2.lags=1, d=1, threshold=0.2)
#' fitrob <- TARMA.fit(x,tar1.lags=1, tar2.lags=1, tma1.lags=1, tma2.lags=1, d=1,
#'             method='robust',alpha=0.7,qu=c(0.1,0.95), threshold=0.2)
#'
#' par.true <- c(0.6,0.6,-1,0.4,-0.7,0.5)
#' pnames   <- c("int.1", "ar1.1", "int.2", "ar2.1", "ma1.1", "ma2.1")
#' names(par.true) <- pnames
#'
#' par.ls  <- round(fitls$fit$coef,2)  # Least Squares
#' par.rob <- round(fitrob$fit$coef,2) # robust
#'
#' rbind(par.true,par.ls,par.rob)

## ***************************************************************************

TARMA.fit <- function(x, tar1.lags = c(1), tar2.lags = c(1), tma1.lags = c(1),
                         tma2.lags = c(1), threshold = NULL, d = 1, pa = 0.25, pb = 0.75,
                         method=c("L-BFGS-B","solnp","lbfgsb3c","robust","trimmed"),
                         alpha=0, qu=c(0.05,0.95), innov = c('norm','student'), optim.control = list(trace=0), 
                         irls.control = list(maxiter=100, tol=1e-4),  ...) {
  
  method   <- match.arg(method)
  innov    <- match.arg(innov)
  n        <- length(x)
  n.tar1   <- length(tar1.lags)
  n.tar2   <- length(tar2.lags)
  t.lags   <- sort(unique(c(tar1.lags,tar2.lags)))
  n.tma1   <- length(tma1.lags)
  n.tma2   <- length(tma2.lags)
  m.lags   <- sort(unique(c(tma1.lags,tma2.lags)))
  if(any(t.lags<=0)||any(m.lags<=0)) stop('lags must be positive integers')
  
  if(!any(is.ts(x)|is.zoo(x))){x <- ts(x)}
  ix <- index(x)
  p1 <- n.tar1+1 # dimension of the TAR vector (1st regime)
  p2 <- n.tar2+1 # dimension of the TAR vector (2nd regime)
  q1 <- n.tma1   # dimension of the TMA vector (1st regime)
  q2 <- n.tma2   # dimension of the TMA vector (2nd regime)
  n.tarp <- p1+p2 # total number of TAR parameters
  n.tmap <- q1+q2 # total number of TMA parameters
  ptot   <- n.tarp+n.tmap # total number of parameters
  th     <- rep(0,ptot)  # parameter vector
  k      <- max(t.lags,m.lags,d,na.rm=TRUE) # number of discarded observations
  neff   <- n-k # actual sample size
  xth    <- c(rep(mean(x),d),x[1:(n-d)])  # threshold variable
  indg   <- (k+1):n                 # default, non robust set of indices used to compute the OLS target function and gradient
  wt     <- rep(1,n)                # default, non robust set of weights used to compute the OLS target function and gradient
  ng     <- length(indg)
  
  ## ************************************************
  ## ** LS criteria and their derivatives
  ## ************************************************
  
  tarmals.R <- function(th){
    # TARMA OLS target function (RECURSIVE VERSION)
    # th = (phi1,phi2,th1,th2) parameter vector
    # indg  : external set of indexes (subset of 1:n) for trimmed estimation
    phi1 <- th[1:p1]
    phi2 <- th[(p1+1):(p1+p2)]
    th1  <- th[(p1+p2+1):(p1+p2+q1)]
    th2  <- th[(p1+p2+q1+1):(p1+p2+q1+q2)]
    eps <- rep(0,n)
    L <- 0
    for(i in indg){
      eps[i] <- x[i] - (phi1%*%c(1,x[i-tar1.lags])+th1%*%eps[i-tma1.lags])*Ir[i]  - (phi2%*%(c(1,x[i-tar2.lags]))+th2%*%eps[i-tma2.lags])*(1-Ir[i])
      L <- L + eps[i]^2
    }
    return(L)
  }
  ## ************************************************
  tarmals.f <- function(th){
    # TARMA OLS target function (Fortran)
    # th = (phi1,phi2,th1,th2) parameter vector
    #   SUBROUTINE tarmaLS(x,nx,th,k,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,L)
    rss <- .Fortran('tarmals',as.double(x),as.integer(n), as.double(th),
                    as.integer(k),as.integer(tar1.lags),as.integer(p1),as.integer(tar2.lags),as.integer(p2),
                    as.integer(tma1.lags),as.integer(q1),as.integer(tma2.lags),as.integer(q2),
                    as.integer(Ir), L=double(1),PACKAGE='tseriesTARMA')$L
    #    as.integer(Ir), L=double(1),PACKAGE="tseriesTARMA")$L
    return(rss)
  }
  ## ************************************************
  
  tarma.eps.R <- function(x,Ir,th,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags){
    # computes the residuals of a TARMA fit
    # th = (phi1,phi2,th1,th2) parameter vector
    n  <- length(x)
    p1 <- length(tar1.lags)+1 # dimension of the TAR vector (1st regime)
    p2 <- length(tar2.lags)+1 # dimension of the TAR vector (2nd regime)
    q1 <- length(tma1.lags)   # dimension of the TMA vector (1st regime)
    q2 <- length(tma2.lags)   # dimension of the TMA vector (2nd regime)
    phi1 <- th[1:p1]
    phi2 <- th[(p1+1):(p1+p2)]
    th1  <- th[(p1+p2+1):(p1+p2+q1)]
    th2  <- th[(p1+p2+q1+1):(p1+p2+q1+q2)]
    eps   <- x
    eps[] <- 0
    for(i in indg){
      eps[i] <- x[i] - (phi1%*%c(1,x[i-tar1.lags])+th1%*%eps[i-tma1.lags])*Ir[i]  - (phi2%*%(c(1,x[i-tar2.lags]))+th2%*%eps[i-tma2.lags])*(1-Ir[i])
    }
    return(eps)
  }
  ## ************************************************
  
  Deps.R <- function(x,eps,Ir,th,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags){
    ## gradient of eps wrt the TARMA parameters 
    n  <- length(x)
    p1 <- length(tar1.lags)+1 # dimension of the TAR vector (1st regime)
    p2 <- length(tar2.lags)+1 # dimension of the TAR vector (2nd regime)
    q1 <- length(tma1.lags)   # dimension of the TMA vector (1st regime)
    q2 <- length(tma2.lags)   # dimension of the TMA vector (2nd regime)
    phi1 <- th[1:p1]
    phi2 <- th[(p1+1):(p1+p2)]
    th1  <- th[(p1+p2+1):(p1+p2+q1)]
    th2  <- th[(p1+p2+q1+1):(p1+p2+q1+q2)]
    ptot <- length(th)
    deps <- matrix(0,n,ptot)
    pnames <- c(paste('phi1',c(0,tar1.lags),sep='.'),paste('phi2',c(0,tar2.lags),sep='.'),
                paste('th1',tma1.lags,sep='.'),paste('th2',tma2.lags,sep='.'))
    colnames(deps) <- pnames
    for(i in indg){
      deps[i,1] <- -Ir[i]*(1 + th1%*%deps[i-tma1.lags,1]) - (1-Ir[i])*(th2%*%deps[i-tma2.lags,1])
      deps[i,2:p1] <- -Ir[i]*(x[i-tar1.lags] + th1%*%deps[i-tma1.lags,2:p1]) -
        (1-Ir[i])*(th2%*%deps[i-tma2.lags,2:p1])
      deps[i,(p1+1)] <- -Ir[i]*(th1%*%deps[i-tma1.lags,(p1+1)]) - (1-Ir[i])*(1+th2%*%deps[i-tma2.lags,(p1+1)])
      deps[i,(p1+2):(p1+p2)] <- -Ir[i]*(th1%*%deps[i-tma1.lags,(p1+2):(p1+p2)]) -
        (1-Ir[i])*(x[i-tar2.lags]+th2%*%deps[i-tma2.lags,(p1+2):(p1+p2)])
      deps[i,(p1+p2+1):(p1+p2+q1)] <- -Ir[i]*(eps[i-tma1.lags]+th1%*%deps[i-tma1.lags,(p1+p2+1):(p1+p2+q1)]) -
        (1-Ir[i])*(th2%*%deps[i-tma2.lags,(p1+p2+1):(p1+p2+q1)])
      deps[i,(p1+p2+q1+1):(p1+p2+q1+q2)] <- -Ir[i]*(th1%*%deps[i-tma1.lags,(p1+p2+q1+1):(p1+p2+q1+q2)]) -
        (1-Ir[i])*(eps[i-tma2.lags]+th2%*%deps[i-tma2.lags,(p1+p2+q1+1):(p1+p2+q1+q2)])
    }
    return(deps)
  }
  ## ************************************************
  
  D2eps.R <- function(x,deps,Ir,th,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags){
    ## Hessian of eps wrt the TARMA parameters 
    vaprodsum <- function(v,A){
      # multiplies the scalar element of the vector v
      # by the sheet of the array A and sums them
      # the result is a matrix with the first two dimensions of A
      if(length(v)==1){
        A   <- A[,,1]
        res <- v*A
      }else{
        dimA <- dim(A)
        k    <- dimA[3]
        res <- matrix(0,dimA[1],dimA[2])
        for(i in 1:k){
          res <- res + v[i]*A[,,i]
        }
      }
      return(res)
    }
    ## Hessian of eps
    n  <- length(x)
    p1 <- length(tar1.lags)+1 # dimension of the TAR vector (1st regime)
    p2 <- length(tar2.lags)+1 # dimension of the TAR vector (2nd regime)
    q1 <- length(tma1.lags)   # dimension of the TMA vector (1st regime)
    q2 <- length(tma2.lags)   # dimension of the TMA vector (2nd regime)
    phi1 <- th[1:p1]
    phi2 <- th[(p1+1):(p1+p2)]
    th1  <- th[(p1+p2+1):(p1+p2+q1)]
    th2  <- th[(p1+p2+q1+1):(p1+p2+q1+q2)]
    
    ptot <- length(th)
    d2eps <- array(0,dim=c(ptot,ptot,n))
    pnames <- c(paste('phi1',c(0,tar1.lags),sep='.'),paste('phi2',c(0,tar2.lags),sep='.'),
                paste('th1',tma1.lags,sep='.'),paste('th2',tma2.lags,sep='.'))
    dimnames(d2eps) <- list(pnames,pnames,1:n)
    for(i in indg){
      
      d2eps[(p1+p2+1):(p1+p2+q1),1:p1,i] <- Ht1p1 <- # dth1 dphi1
        -Ir[i]*(deps[i-tma1.lags,1:p1] + vaprodsum(th1,d2eps[(p1+p2+1):(p1+p2+q1),1:p1,i-tma1.lags,drop=FALSE])) -
        (1-Ir[i])*(vaprodsum(th2,d2eps[(p1+p2+1):(p1+p2+q1),1:p1,i-tma2.lags,drop=FALSE]))  
      
      d2eps[1:p1,(p1+p2+1):(p1+p2+q1),i] <- t(Ht1p1) # dphi1 dth1
      
      d2eps[(p1+p2+1):(p1+p2+q1),(p1+1):(p1+p2),i] <- Ht1p2 <-
        -Ir[i]*(deps[i-tma1.lags,(p1+1):(p1+p2)] + vaprodsum(th1,d2eps[(p1+p2+1):(p1+p2+q1),(p1+1):(p1+p2),i-tma1.lags,drop=FALSE])) -
        (1-Ir[i])*(vaprodsum(th2,d2eps[(p1+p2+1):(p1+p2+q1),(p1+1):(p1+p2),i-tma2.lags,drop=FALSE]))  # dth1 dphi2
      
      d2eps[(p1+1):(p1+p2),(p1+p2+1):(p1+p2+q1),i] <- t(Ht1p2) # dphi2 dth1
      
      d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),1:p1,i] <- Ht2p1 <- 
        -Ir[i]*(vaprodsum(th1,d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),1:p1,i-tma1.lags,drop=FALSE])) -
        (1-Ir[i])*(deps[i-tma2.lags,1:p1]+ vaprodsum(th2,d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),1:p1,i-tma2.lags,drop=FALSE]))  # dth2 dphi1
      
      d2eps[1:p1,(p1+p2+q1+1):(p1+p2+q1+q2),i] <- t(Ht2p1)
      
      d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),(p1+1):(p1+p2),i] <- Ht2p2 <- 
        -Ir[i]*(vaprodsum(th1,d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),(p1+1):(p1+p2),i-tma1.lags,drop=FALSE])) -
        (1-Ir[i])*(deps[i-tma2.lags,(p1+1):(p1+p2)]+
                     vaprodsum(th2,d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),(p1+1):(p1+p2),i-tma2.lags,drop=FALSE]))  # dth2 dphi2
      
      d2eps[(p1+1):(p1+p2),(p1+p2+q1+1):(p1+p2+q1+q2),i] <- t(Ht2p2)  # dphi2 dth2
      
      d2eps[(p1+p2+1):(p1+p2+q1),(p1+p2+1):(p1+p2+q1),i] <- -Ir[i]*(
        deps[i-tma1.lags,(p1+p2+1):(p1+p2+q1)] + t(deps[i-tma1.lags,(p1+p2+1):(p1+p2+q1)])
        + vaprodsum(th1,d2eps[(p1+p2+1):(p1+p2+q1),(p1+p2+1):(p1+p2+q1),i-tma1.lags,drop=FALSE])) 
      - (1-Ir[i])*(vaprodsum(th2,d2eps[(p1+p2+1):(p1+p2+q1),(p1+p2+1):(p1+p2+q1),i-tma2.lags,drop=FALSE])) # dth1 dth1
      
       d2eps[(p1+p2+1):(p1+p2+q1),(p1+p2+q1+1):(p1+p2+q1+q2),i] <- Ht1t2 <- -Ir[i]*(
         deps[i-tma1.lags,(p1+p2+q1+1):(p1+p2+q1+q2)]
         + vaprodsum(th1,d2eps[(p1+p2+1):(p1+p2+q1),(p1+p2+q1+1):(p1+p2+q1+q2),i-tma1.lags,drop=FALSE]))
       - (1-Ir[i])*(t(deps[i-tma2.lags,(p1+p2+1):(p1+p2+q1)]) +
                      vaprodsum(th2,d2eps[(p1+p2+1):(p1+p2+q1),(p1+p2+q1+1):(p1+p2+q1+q2),i-tma2.lags,drop=FALSE])) # dth1 dth2
       
       d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),(p1+p2+1):(p1+p2+q1),i] <- t(Ht1t2) # dth2 dth1
      
      d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),(p1+p2+q1+1):(p1+p2+q1+q2),i] <- -Ir[i]*(
        vaprodsum(th1,d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),(p1+p2+q1+1):(p1+p2+q1+q2),i-tma1.lags,drop=FALSE]))
      - (1-Ir[i])*(deps[i-tma2.lags,(p1+p2+q1+1):(p1+p2+q1+q2)] + t(deps[i-tma2.lags,(p1+p2+q1+1):(p1+p2+q1+q2)])
                   + vaprodsum(th2,d2eps[(p1+p2+q1+1):(p1+p2+q1+q2),(p1+p2+q1+1):(p1+p2+q1+q2),i-tma2.lags,drop=FALSE])) # dth2 dth2
    }
    return(d2eps)
  }
  ## ************************************************
  #DLeps.R <- function(th){
  ### Analytical Gradient of the TARMA least squares criterion
  #    eps  <- as.double(tarma.eps.R(x,Ir,th,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags))
  #    deps <- Deps.R(x,eps,Ir,th,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags)
  #    res  <- colSums(2*eps*deps)
  #    return(res)
  #}
  ## ************************************************
  DLeps.f <- function(th){
    ## Analytical Gradient of the TARMA least squares criterion
    ## Fortran version
    #SUBROUTINE tarmaDLS(x,nx,th,k,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,DL)
    ptot <- length(th)
    res <- .Fortran('tarmadls',as.double(x),as.integer(n), as.double(th),
                    as.integer(k),as.integer(tar1.lags),as.integer(p1),as.integer(tar2.lags),as.integer(p2),
                    as.integer(tma1.lags),as.integer(q1),as.integer(tma2.lags),as.integer(q2),
                    as.integer(Ir), DL=double(ptot),PACKAGE='tseriesTARMA')$DL
    return(res)
  }
  ## ************************************************
  
  if(method=='robust'| method=='trimmed'){
    ## ************************************************
    #    tarmalsw.R <- function(th,wt){
    #    # TARMA OLS target function (ROBUST VERSION)
    #    # th    :  (phi1,phi2,th1,th2) parameter vector
    #    # wt    : weights (that depend upon alpha)
    #    # indg  : external set of indexes (subset of 1:n) for trimmed estimation
    #        phi1 <- th[1:p1]
    #        phi2 <- th[(p1+1):(p1+p2)]
    #        th1  <- th[(p1+p2+1):(p1+p2+q1)]
    #        th2  <- th[(p1+p2+q1+1):(p1+p2+q1+q2)]
    #        eps <- rep(0,n)
    #        L <- 0
    #        for(i in indg){
    #            eps[i] <- x[i] - (phi1%*%c(1,x[i-tar1.lags])+th1%*%eps[i-tma1.lags])*Ir[i]  - (phi2%*%(c(1,x[i-tar2.lags]))+th2%*%eps[i-tma2.lags])*(1-Ir[i])
    #            L <- L + wt[i]*eps[i]^2
    #        }
    #        return(L)
    #    }
    ## ************************************************
    tarmalsw.f <- function(th,wt){
      # TARMA OLS target function (ROBUST VERSION, Fortran)
      # wt    : weights (that depend upon alpha)
      # indg  : external set of indexes (subset of 1:n) for trimmed estimation
      
      # SUBROUTINE tarmaLSW(x,nx,th,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,wt,indg,ng,L)
      rss <- .Fortran('tarmalsw',as.double(x),as.integer(n), as.double(th),
                      as.integer(tar1.lags),as.integer(p1),as.integer(tar2.lags),as.integer(p2),
                      as.integer(tma1.lags),as.integer(q1),as.integer(tma2.lags),as.integer(q2),
                      as.integer(Ir), as.double(wt), as.integer(indg),as.integer(ng),L=double(1),PACKAGE='tseriesTARMA')$L
      #    ,PACKAGE="tseriesTARMA")$L
      return(rss)
    }
    ## ************************************************
    #    DLepsw.R <- function(th,wt){
    #    ## Gradient of the robust TARMA least squares criterion tarmalsw.R
    #    # wt    : weights (that depend upon alpha)
    #        eps  <- as.double(tarma.eps.R(x,Ir,th,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags))
    #        deps <- Deps.R(x,eps,Ir,th,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags)
    #        res  <- colSums(2*wt*eps*deps)
    #        return(res)
    #    }
    ## ************************************************
    DLepsw.f <- function(th,wt){
      ## Analytical Gradient of the robust TARMA least squares criterion
      ## Fortran version
      # wt    : weights (that depend upon alpha)
      # SUBROUTINE tarmaDLSW(x,nx,th,tlag1,p1,tlag2,p2,mlag1,q1,mlag2,q2,Ir,wt,indg,ng,DL)
      ptot <- length(th)
      res <- .Fortran('tarmadlsw',as.double(x),as.integer(n), as.double(th),
                      as.integer(tar1.lags),as.integer(p1),as.integer(tar2.lags),as.integer(p2),
                      as.integer(tma1.lags),as.integer(q1),as.integer(tma2.lags),as.integer(q2),
                      as.integer(Ir), as.double(wt), as.integer(indg),as.integer(ng), DL=double(ptot),PACKAGE='tseriesTARMA')$DL
      return(res)
    }
    ## ****************************************************
    ## Gradient of the rho function
    ## ****************************************************
    Drho <- function(type,alpha){
      # first derivative of rho wrt to eps for different densities for the innovations
      if(type=='norm'){
        function(eps,alpha){
          s2t <- mean(eps^2)
#          (2*pi)^(-alpha/2)*s2t^(-alpha/2-1)*c(exp(-alpha*eps^2/(2*s2t)))*eps
          (2*pi)^(-alpha/2)*s2t^(-1)*c(exp(-alpha*eps^2/(2*s2t)))*eps
        }
      # }else if(type=='sn'){
      #   if(alpha>0){
      #     function(eps,alpha){
      #       fit  <- selm(eps ~ 1, family="SN")
      #       pa   <- fit@param$dp # SN parameters (xi,omega,a)
      #       #         pa[1] <- -pa[2]*sqrt(2/pi)*pa[3]/sqrt(1+pa[3]^2) # to get zero mean
      #       ec    <- (eps-pa[1])/pa[2]    
      #       #      drho <- pa[2]^(-3) * exp(-(1+pa[3]^2)/2*ec^2) * (2*pi)^(-alpha/2) * (exp(-ec^2/2) * 2*(pnorm(pa[3]*ec))/pa[2])^(alpha-1)*
      #       #        (-pa[3]*sqrt(2/pi)*pa[2] + exp(pa[3]^2/2 * ec^2) * (eps-pa[1])* 2*pnorm(pa[3]*ec))
      #       drho <- pa[2]^(-3) * (2*pi)^(-alpha/2)* (2*pnorm(pa[3]*ec)/pa[2])^(alpha) * exp(-alpha*ec^2/2)  *
      #         (pa[2]*(eps-pa[1]) - pa[3]*sqrt(2/pi)*pa[2]^2 *exp(-pa[3]^2 * ec^2/2)*(2*pnorm(pa[3]*ec))^(-1))
      #       return(drho)
      #     }
      #   }else{
      #     function(eps,alpha){
      #       fit  <- selm(eps ~ 1, family="SN")
      #       pa   <- fit@param$dp # SN parameters (xi,omega,a)
      #       ec    <- (eps-pa[1])/pa[2]    
      #       drho <- (eps - pa[1] - (pa[3]* exp(-(pa[3]*ec)^2/2)*sqrt(2/pi)*pa[2])/((2*pnorm(pa[3]*ec))))/pa[2]^2
      #       return(drho)
      #     }
      #   }
      }else if(type=='student'){
        function(eps,alpha){
          vx  <- var(eps)
          sdf <- ifelse(vx<=1.3,20,2/(1-1/vx))
          f0  <- tryCatch(fitdist(eps, distr='t', method = "mge", start=list(df=sdf)),error=function(x){list(estimate=sdf)})
          sdf2 <- as.double(f0$estimate)
          f1   <- tryCatch(fitdist(eps, distr='t', method = "mle",start=list(df=sdf2)),error=function(x){list(estimate=sdf2)})
          nu  <- f1$estimate
          gamma((nu+1)/2)^alpha * gamma(nu/2)^(-alpha)*(nu*pi)^(-alpha/2)*(1+eps^2/nu)^(-alpha*(nu+1)/2)*(nu+1)/(2*(nu+eps^2))*2*eps
        }
      }
    }
    
    ## ********************************************************
    D2rho <- function(type,alpha){
      # second derivative of rho wrt to eps for different densities for the innovations
      if(type=='norm'){
        function(eps,alpha){
          s2t <- mean(eps^2)
#          (2*pi)^(-alpha/2)*s2t^(-alpha/2-2)* c(exp(-alpha*eps^2/(2*s2t)))*(-alpha*eps^2 + s2t)
          (2*pi)^(-alpha/2)*s2t^(-2)* c(exp(-alpha*eps^2/(2*s2t)))*(-alpha*eps^2 + s2t)
        }
      # }else if(type=='sn'){
      #   if(alpha>0){
      #     function(eps,alpha){
      #       fit  <- selm(eps ~ 1, family="SN")
      #       pa   <- fit@param$dp # SN parameters (xi,omega,a)
      #       #          pa[1] <- -pa[2]*sqrt(2/pi)*pa[3]/sqrt(1+pa[3]^2) # to get zero mean
      #       ec    <- (eps-pa[1])/pa[2]
      #       Phi   <-  pnorm(pa[3]*ec)
      #       out <- ((pa[3]^2*(1-alpha))/(2*(pa[2]*Phi)^2) * 2^(-alpha/2) *exp(-(pa[3]*ec)^2)*pi^(-1-alpha/2) +
      #                 pa[3]*(2^((1-alpha)/2)*pi^(-(1+alpha)/2))/(pa[2]^2*2*Phi) * ec*(pa[3]^2 + 2*alpha)* exp(-(pa[3]*ec)^2/2)
      #               + pa[2]^(-2)*(2*pi)^(-alpha/2)*(1-alpha*ec^2)) *  (exp(-ec^2/2) *2*Phi/pa[2])^(alpha)
      #       return(out)
      #     }
      #   }else{
      #     function(eps,alpha){
      #       fit  <- selm(eps ~ 1, family="SN")
      #       pa   <- fit@param$dp # SN parameters (xi,omega,a)
      #       ec   <- (eps-pa[1])/pa[2]    
      #       out <- pa[2]^(-2)*(1+(1*pa[3]^2*exp(-(pa[3]*ec)^2))/(2*pi*pnorm(pa[3]*ec)^2) + 
      #                            (pa[3]^3*exp(-(pa[3]*ec)^2/2)*sqrt(2/pi)*(eps-pa[1]))/(2*pa[2]*pnorm(pa[3]*ec)))
      #       return(out)
      #     }
      #   }
      }else if(type=='student'){
        function(eps,alpha){
          vx  <- var(eps)
          sdf <- ifelse(vx<=1.3,20,2/(1-1/vx))
          f0  <- tryCatch(fitdist(eps, distr='t', method = "mge", start=list(df=sdf)),error=function(x){list(estimate=sdf)})
          sdf2 <- as.double(f0$estimate)
          f1   <- tryCatch(fitdist(eps, distr='t', method = "mle",start=list(df=sdf2)),error=function(x){list(estimate=sdf2)})
          nu  <- f1$estimate
          -gamma((nu+1)/2)^alpha * gamma(nu/2)^(-alpha)*(nu*pi)^(-alpha/2)*(1+eps^2/nu)^(-alpha*(nu+1)/2)*(nu+1)/(2*(nu+eps^2))*2*
            (-nu+eps^2+(1+nu)*eps^2 * alpha)/((nu+eps^2))
        }
      }
    }
    drho   <- Drho(type=innov,alpha)
    d2rho  <- D2rho(type=innov,alpha)
    ## ***********************************************************
    wtfun <- function(eps,alpha){
      ## function to compute the robustness weights
      out    <- rep(0,length(eps))
      ind    <- which(abs(eps)>0)
      dreps  <- drho(eps,alpha) # first derivative of the rho function w.r.t. eps
      out[ind] <- dreps[ind]/(2*eps[ind])
      #     out      <- ifelse(abs(out)>1,1,out)
      return(out)    
    }
    ## ***********************************************************
    
    irls <-  function(par, optf, par.opt, control=list()){
      # Iteratively Re-weighted Least Squares step
      con <- list(maxiter=100, tol=1e-04)
      con[names(control)] <- control
      maxiter  <- con$maxiter
      tol      <- con$tol
      p.old    <- par
      crit.old <- Inf # starting value for the reduction criterion 
      eps <- tarma.eps.R(x,Ir,par,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags)
      wt  <- wtfun(eps,alpha)
      par.opt$par <- par
      par.opt$wt  <- wt
      for(j in 1:maxiter){
        opt.new      <- tryCatch(do.call(optf,par.opt),
                                 error=function(x){list(par=rep(Inf,ptot),value=Inf, convergence=-1)})
        ##            opt.new  <- do.call(optf,par.opt)
        p.new    <- opt.new$par
        crit.new <- norm((p.new-p.old)/p.old,'2')
        #            cat(crit.new,opt.new$value,'\n')
        ##            if((crit.new > crit.old)|(crit.new < tol)) break
        if(crit.new < tol) break
        if((opt.new$convergence!=0)|(crit.new > crit.old*100)){
          p.new <- p.old + rnorm(ptot,sd=0.05)
        }else{
          opt <- opt.new
          crit.old <- crit.new
        }
        eps <- tarma.eps.R(x,Ir,th=p.new,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags)
        wt  <- wtfun(eps,alpha)
        par.opt$par <- p.new
        par.opt$wt  <- wt
        p.old       <- p.new
      }
      s2t <- sum(eps^2)/(n - k)
      opt$iterations <- j
      opt$sigma2     <- s2t/(1-alpha/2)
      opt$wt         <- wt
      #        cat('IRLS iterations:',j,'------------ \n')
      return(opt)
    }
  }
  
  ## **************************************************************************
  ## Constraints for ergodicity and invertibility
  ## **************************************************************************
  
  if(all(tar1.lags==1)&all(tar2.lags==1)&all(tma1.lags==1)&all(tma2.lags==1)){
    # TARMA(1,1)
    tarma.parcheck <- function(th){
      # constraints for geometric ergodicity and invertibility
      # for a TARMA(1,1) model
      phi1 <- th[2]
      phi2 <- th[4]
      res  <- (phi1<1)&(phi2<1)&((phi1*phi2)<1)&(abs(th[5])<1)&(abs(th[6])<1)
      return(unname(res))
    }
    ifun <- function(th){
      # inequality constraints for geometric ergodicity and invertibility
      # for a TARMA(1,1) model
      phi1 <- th[2]
      phi2 <- th[4]
      res  <- c(phi1*phi2)
    }
    ilow.b <- -Inf
    iupp.b <-  1
    lowb <- c(-Inf,-Inf,-Inf,-Inf,-1,-1)
    uppb <- c(+Inf,1,+Inf,1,1,1)
    
  }else{
    # generic TARMA
    tarma.parcheck <- function(th){
      # constraints for geometric ergodicity
      # for a TARMA model
      phi1 <- th[2:p1]           # without intercepts
      phi2 <- th[(p1+2):(p1+p2)] # without intercepts
      res  <- (sum(abs(phi1))<1)&(sum(abs(phi2))<1)
      return(unname(res))
    }
    ifun <- function(th){
      # inequality constraints for geometric ergodicity and invertibility
      # for a TARMA model
      phi1 <- th[2:p1]           # without intercepts
      phi2 <- th[(p1+2):(p1+p2)] # without intercepts
      th1  <- th[(p1+p2+1):(p1+p2+q1)]
      th2  <- th[(p1+p2+q1+1):(p1+p2+q1+q2)]
      res  <- c(sum(abs(phi1)),sum(abs(phi2)),sum(abs(th1)),sum(abs(th2)))
    }
    ilow.b <- rep(0,4)
    iupp.b <- rep(1,4)
    # parameter bounds
    lowb      <- c(-Inf,rep(-1,p1-1),-Inf,rep(-1,p2-1),rep(-1,q1),rep(-1,q2))    # TARMA(p,q)
    uppb      <- c(+Inf,rep( 1,p1-1),+Inf,rep(+1,p2-1),rep(+1,q1),rep(+1,q2))
  }
  ## *****************************************************
  ## OPTIMIZATION METHODS ********************************
  ## *****************************************************
  
  if(method=='L-BFGS-B'){
    optfun <- 'optim'
    fobj   <- tarmals.f
    gobj   <- DLeps.f
    param  <- list(fn=fobj,gr=gobj,method="L-BFGS-B", lower=lowb,upper=uppb,
                   control=optim.control)
  }else if(method=='solnp'){
    optfun <- 'solnp'
    fobj   <- tarmals.f
    gobj   <- DLeps.f
    param  <- list(fun=fobj, eqfun = NULL, eqB = NULL, ineqfun = ifun, ineqLB = ilow.b,
                   ineqUB = iupp.b, LB = lowb, UB = uppb, control=optim.control)
  }else if(method=='robust'){
    optfun <- 'irls'
    optf   <- 'optim'
    fobj   <- tarmalsw.f
    gobj   <- DLepsw.f
    par.opt  <- list(fn=fobj,gr=gobj,method="L-BFGS-B", lower=lowb,upper=uppb,
                     control=optim.control)
    param <- list(optf=optf,par.opt=par.opt, control=irls.control)
    # ------------
  }else if(method=='lbfgsb3c'){
    optfun <- 'lbfgsb3c'
    fobj   <- tarmals.f
    gobj   <- DLeps.f
    param  <- list(fn=fobj,gr=gobj,lower=lowb,upper=uppb,control=optim.control)
    #}else if(method=='robustc'){
    #    optfun <- 'lbfgsb3c'
    #    fobj   <- tarmalsw.f
    #    gobj   <- DLepsw.f
    #    param  <- list(fn=fobj,gr=gobj, lower=lowb,upper=uppb,
    #     control=optim.control)
    #}else if(method=='rsolnp'){
    #    optfun <- 'solnp'
    #    fobj   <- tarmalsw.f
    #    gobj   <- DLepsw.f
    #    param  <- list(fun=fobj, eqfun = NULL, eqB = NULL, ineqfun = ifun, ineqLB = ilow.b,
    #           ineqUB = iupp.b, LB = lowb, UB = uppb,control=optim.control)
  }else if(method=='trimmed'){
    quant  <- quantile(x,probs=qu)
    indg   <- intersect(indg,which((quant[1]<x)&(x<quant[2]))) # trimmed set of indices 
    optfun <- 'optim'
    fobj   <- tarmalsw.f
    gobj   <- DLepsw.f
    param  <- list(fn=fobj,gr=gobj,method="L-BFGS-B", lower=lowb,upper=uppb,
                   control=optim.control)
  }
  ## *********************************************
  if(is.null(threshold)){ # estimates the threshold
    a         <- ceiling((neff-1)*pa)
    b         <- floor((neff-1)*pb)
    thd.v     <- sort(xth)
    thd.range <- thd.v[a:b]
    nr        <- length(thd.range)
    pars      <- matrix(NA,nr,ptot)
    rss.g     <- Inf
    par.g     <- rep(NA,ptot)
    rssv      <- rep(NA,nr)
    ind       <- trunc(quantile(1:nr,prob=seq(0,1,by=0.1)))
    res       <- matrix(NA,nr,ptot+1)
    x.t      <- x       # for the initial fit (either trimmed or not)   
    if((method=='robust')&(alpha>0)){
      quant    <- quantile(x,probs=c(qu[1],qu[2]))
      indg     <- intersect(indg,which((quant[1]<x)&(x<quant[2]))) # trimmed set of indices for robust initial fit
      ng       <- length(indg)
      indb     <- setdiff(1:n,indg)
      x.t[indb]<- NA; #median(x[indg])
    }
    fit0      <- TARMA.fit2(x.t, ma.ord = max(m.lags), ar.lags = NULL,
                            tar1.lags = tar1.lags, tar2.lags=tar2.lags, threshold=quantile(x.t,pa,na.rm=TRUE)
                            , d = d, include.int = TRUE)
    ar1      <- fit0$phi1[c(1,tar1.lags+1)]
    ar2      <- fit0$phi2[c(1,tar2.lags+1)]
    ma1      <- fit0$theta1[tma1.lags]
    ma2      <- fit0$theta2[tma2.lags]
    theta0   <- c(ar1,ar2,ma1,ma2)+rnorm(ptot,sd=0.01)
    thd      <- thd.range[1]
    Ir       <- (xth<=thd)
    sigma2   <- fit0$fit$sigma2
    if((method=='robust')&(alpha>0)){ # initial robust fit using trimmed LS
      optfun.s <- 'optim'
      fobj.s   <- tarmalsw.f
      gobj.s   <- DLepsw.f
      param.s  <- list(par=theta0, fn=fobj.s,gr=gobj.s,method="L-BFGS-B",lower=lowb,upper=uppb,
                       control=optim.control,wt=wt)
      opt      <- tryCatch(do.call(optfun.s,param.s),
                           error=function(x){list(par=rep(NA,ptot),value=Inf, convergence=-1)})
      theta0   <- opt$par
      indg   <- (k+1):n  # 
    }
    param$par <- theta0
    opt       <- tryCatch(do.call(optfun,param),
                          error=function(x){list(par=rep(NA,ptot),value=Inf, convergence=-1)})
    ntry <- 0
    #        while(opt$convergence!=0){
    while(is.null(opt$convergence)||opt$convergence!=0){ # workaround for the bug in lbfgsb3c
      ntry      <- ntry+1
      theta0    <- theta0+rnorm(ptot,sd=0.05)
      param$par <- theta0
      opt       <- tryCatch(do.call(optfun,param),
                            error=function(x){list(par=rep(NA,ptot),value=Inf, convergence=-1)})
      if(ntry==50) stop('Initial fit failed')
    }
    res[1,1]   <- opt$value[1]
    res[1,-1]  <- opt$par
    theta0     <- opt$par
    for(i in 2:nr){
      thd       <- thd.range[i]
      Ir        <- (xth<=thd)
      param$par <- theta0
      opt       <- tryCatch(do.call(optfun,param),
                            error=function(x){list(par=rep(NA,ptot),value=Inf, convergence=-1)})
      #            if(opt$convergence!=0){
      ntry <- 0
      while(is.null(opt$convergence)||opt$convergence!=0){ # workaround for the bug in lbfgslb3c
        if(ntry>5){theta0 <- res[i-2,-1]; break}
        #              theta0    <- theta0 + rnorm(ptot,sd=0.05)
        theta0    <- res[i-2,-1] + rnorm(ptot,sd=0.03) # uses the last good parameter
        param$par <- theta0
        opt       <- tryCatch(do.call(optfun,param),
                              error=function(x){list(par=rep(NA,ptot),value=Inf, convergence=-1)})
        ntry <- ntry + 1
#        cat('ntry: ',ntry,'\n')
      }
      theta0 <- opt$par
      res[i,1]  <- opt$value[1]
      res[i,-1] <- opt$par
    }
    igood <- which.min(res[,1])
    thd   <- thd.range[igood]
    rssv  <- res[,1]
    rss   <- res[igood,1]
    par.g <- res[igood,-1]
    Ir    <- (xth<= thd)
  }else{ # fixed threshold -----------------------------------------------
    thd       <- threshold
    Ir        <- (xth<= thd)
    rssv      <- NA
    thd.range <- NA
    if((method=='robust')&(alpha>0)){ # initial robust fit using trimmed LS
      quant    <- quantile(x,probs=c(qu[1],qu[2]))
      indg     <- intersect(indg,which((quant[1]<x)&(x<quant[2]))) # trimmed set of indices for robust initial fit
      ng       <- length(indg)
      indb     <- setdiff(1:n,indg)
      x.t      <- x                 # trimmed x
      x.t[indb]<- NA #median(x[indg])
      #          xth.t    <- c(rep(median(x.t),d),x.t[1:(n-d)])  # trimmed threshold variable
      fit0     <- TARMA.fit2(x.t, ma.ord=max(m.lags), ar.lags=NULL,
                             tar1.lags=tar1.lags, tar2.lags=tar2.lags, threshold=thd, d=d,include.int=TRUE)
      ar1       <- fit0$phi1[c(1,tar1.lags+1)]
      ar2       <- fit0$phi2[c(1,tar2.lags+1)]
      ma1       <- fit0$theta1[tma1.lags]
      ma2       <- fit0$theta2[tma2.lags]
      sigma2    <- fit0$fit$sigma2
      theta0    <- c(ar1,ar2,ma1,ma2)+rnorm(ptot,sd=0.05)
      optfun.s <- 'optim'
      fobj.s   <- tarmalsw.f
      gobj.s   <- DLepsw.f
      param.s  <- list(par=theta0, fn=fobj.s,gr=gobj.s,method="L-BFGS-B", lower=lowb,upper=uppb,
                       control=optim.control,wt=wt)
      opt      <- tryCatch(do.call(optfun.s,param.s),
                           error=function(x){list(par=rep(NA,ptot),value=Inf, convergence=-1)})
      theta0   <- opt$par
      indg     <- (k+1):n  # 
     }else{
      fit0      <- TARMA.fit2(x, ma.ord=max(m.lags), ar.lags=NULL,
                              tar1.lags=tar1.lags,tar2.lags=tar2.lags,threshold=thd, d=d,include.int=TRUE)
      ar1       <- fit0$phi1[c(1,tar1.lags+1)]
      ar2       <- fit0$phi2[c(1,tar2.lags+1)]
      ma1       <- fit0$theta1[tma1.lags]
      ma2       <- fit0$theta2[tma2.lags]
      sigma2    <- fit0$fit$sigma2
      theta0    <- c(ar1,ar2,ma1,ma2)+rnorm(ptot,sd=0.05)
    }
    param$par <- theta0
    opt       <- do.call(optfun,param)
    rss       <- opt$value[1]
    par.g     <- opt$par
  }
  ## Standard Errors *******************************************************
  innovpar <- NULL # parameter vector for the innovation density (NULL by default)
  wt    <- NULL    # robust weights   
  eps   <- tarma.eps.R(x,Ir,th=par.g,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags)
  deps  <- Deps.R(x,eps,Ir,th=par.g,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags)   # gradient of eps wrt the TARMA parameters
  d2eps <- D2eps.R(x,deps,Ir,th=par.g,k,tar1.lags,tar2.lags,tma1.lags,tma2.lags) # Hessian of eps wrt the TARMA parameters
  eps   <- c(eps[-(1:k)])
  Ir    <- Ir[-(1:k)]
  deps  <- unclass(deps[-(1:k),])
  d2eps <- unclass(d2eps[,,-(1:k)])
  if((method=='robust')|(method=='robustc')|(method=='rsolnp')){
    wt <- opt$wt
    ## Sandwich estimator of the variance/covariance matrix of the estimates
    dreps  <- unclass(drho(eps,alpha))  # first derivative of the rho function w.r.t. eps
    d2reps <- unclass(d2rho(eps,alpha)) # second derivative of the rho function w.r.t. eps
    H   <- matrix(0,ptot,ptot) # Second derivatives 
    J   <- matrix(0,ptot,ptot) # product of the first derivatives 
    for(i in 1:(n-k)){
      deps2i <- deps[i,]%*%t(deps[i,])
      J   <- J  + dreps[i]^2*deps2i
      H   <- H + d2reps[i]*deps2i + dreps[i]*d2eps[,,i]
    }
    if(innov=='student'){
      vx  <- var(eps)
      sdf <- ifelse(vx<=1.3,20,2/(1-1/vx))
      f0  <- tryCatch(fitdist(eps, distr='t', method = "mge", start=list(df=sdf)),error=function(x){list(estimate=sdf)})
      sdf2 <- as.double(f0$estimate)
      f1   <- tryCatch(fitdist(eps, distr='t', method = "mle",start=list(df=sdf2)),error=function(x){list(estimate=sdf2)})
      innovpar <- f1$estimate
      names(innovpar) <- 'df'
    }
  }else{
    D1L   <- 2*c(eps)*deps
    H     <- 2*t(deps)%*%deps + apply(2*array(rep(eps,each=ptot^2),dim=c(ptot,ptot,(n-k)))*d2eps,FUN=sum,MARGIN=c(1,2))
    J     <- t(D1L)%*%D1L
  }
  J   <- J/(n-k)
  H   <- H/(n-k)
  Hi  <- solve(forceSymmetric(H))
  sigmai  <- Hi%*%J%*%Hi
  s2hat   <- mean(eps^2)
  eps     <- as.ts(zoo(eps,order.by = ix[-(1:k)]))   
  n1      <- sum(Ir)
  n2      <- sum(1-Ir)
  np1     <- p1+q1
  np2     <- p2+q2
  s2hat.1 <- sum(Ir*eps^2)/n1
  s2hat.2 <- sum((1-Ir)*eps^2)/n2
  aic     <- n1*log(s2hat.1) + n2*log(s2hat.2) + neff*(1+log(2*pi)) + 2*(np1+np2+1)
  bic     <- n1*log(s2hat.1) + n2*log(s2hat.2) + neff*(1+log(2*pi)) + log(neff)*(np1+np2+1)
  vcoef   <- sigmai/(n-k)
  se      <- sqrt(diag(vcoef))
  ntar1   <- c('int.1',paste('ar1',tar1.lags,sep='.'))
  ntar2   <- c('int.2',paste('ar2',tar2.lags,sep='.'))
  nma1    <- paste('ma1',tma1.lags,sep='.')
  nma2    <- paste('ma2',tma2.lags,sep='.')
  names(par.g) <- names(se) <- c(ntar1,ntar2,nma1,nma2)
  colnames(vcoef) <- rownames(vcoef) <- names(par.g)
  phi1    <- par.g[1:p1]
  phi2    <- par.g[(p1+1):(p1+p2)]
  th1     <- par.g[(p1+p2+1):(p1+p2+q1)]
  th2     <- par.g[(p1+p2+q1+1):(p1+p2+q1+q2)]
  Phi1                   <- rep(0,max(tar1.lags)+1)
  Phi1[c(1,tar1.lags+1)] <- phi1
  names(Phi1)            <- 0:max(tar1.lags)
  Phi2                   <- rep(0,max(tar2.lags)+1)
  Phi2[c(1,tar2.lags+1)] <- phi2
  names(Phi2)            <- 0:max(tar2.lags)
  Th1                    <- rep(0,max(tma1.lags))
  Th1[tma1.lags]         <- th1
  names(Th1)             <- 1:max(tma1.lags)
  Th2                    <- rep(0,max(tma2.lags))
  Th2[tma2.lags]         <- th2
  names(Th2)             <- 1:max(tma2.lags)
  res <- list(fit=list(coef=par.g,sigma2=s2hat,var.coef=vcoef,
                       residuals=eps,nobs=n-k),se=se, thd=thd, aic=aic, bic=bic, rss=rss, rss.v=rssv, thd.range=thd.range,
              d=d, phi1=Phi1, phi2=Phi2, theta1=Th1, theta2=Th2, tlag1=tar1.lags,
              tlag2=tar2.lags, mlag1=tma1.lags, mlag2=tma2.lags, method=method, innov=innov, alpha=alpha, qu=qu,
              call = match.call(), convergence=opt$convergence, innovpar=innovpar, wt=wt)
  if(!tarma.parcheck(par.g)) warning("The parameters are not in the ergodicity region.
             Refit with constrained optimization (method='solnp').")
  if(opt$convergence!=0) warning(paste('convergence code:',opt$convergence))
  class(res) <- 'TARMA'
  return(res)
}

## ****************************************************************************
