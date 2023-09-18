#'  TARMA Modelling of Time Series
#'
#' @description \loadmathjax
#' Maximum Likelihood fit of a two-regime \code{TARMA(p1,p2,q,q)}
#' model with common MA parameters, possible common AR parameters and possible covariates.
#'
#' @param x        A univariate time series.
#' @param ar.lags  Vector of common AR lags. Defaults to \code{NULL}. It can be a subset of lags.
#' @param tar1.lags Vector of AR lags for the lower regime. It can be a subset of \code{1 ... p1 = max(tar1.lags)}.
#' @param tar2.lags Vector of AR lags for the upper regime. It can be a subset of \code{1 ... p2 = max(tar2.lags)}.
#' @param ma.ord   Order of the MA part (also called \code{q} below).
#' @param sma.ord  Order of the seasonal MA part (also called \code{Q} below).
#' @param period   Period of the seasonal MA part (also called \code{s} below).
#' @param estimate.thd Logical. If \code{TRUE} estimates the threshold over the threshold range specified by \code{pa} and \code{pb}.
#' @param threshold Threshold parameter. Used only if \code{estimate.thd = FALSE}.
#' @param d  Delay parameter. Defaults to \code{1}.
#' @param pa  Real number in \code{[0,1]}. Sets the lower limit for the threshold search to the \code{100*pa}-th sample percentile.
#' The default is \code{0.25}
#' @param pb  Real number in \code{[0,1]}. Sets the upper limit for the threshold search to the \code{100*pb}-th sample percentile.
#' The default is \code{0.75}
#' @param thd.var Optional exogenous threshold variable. If \code{NULL} it is set to \code{lag(x,-d)}.
#' If not \code{NULL} it has to be a \code{ts} object.
#' @param include.int Logical. If \code{TRUE} includes the intercept terms in both regimes, a common intercept is included otherwise.
#' @param x.reg  Covariates to be included in the model. These are passed to \code{arima}. If they are not \code{ts} objects they must have the same length as \code{x}.
#' @param optim.control List of control parameters for the optimization method.
#' @param \dots  Additional arguments passed to \code{arima}.
#'
#' @details
#' Fits the following two-regime \code{TARMA} process with optional components: linear \code{AR} part, seasonal \code{MA} and covariates. \cr
#' \mjdeqn{X_{t} = \phi_{0} + \sum_{h \in I} \phi_{h} X_{t-h} + \sum_{l=1}^Q \Theta_{l} \epsilon_{t-ls} + \sum_{j=1}^q \theta_{j} \epsilon_{t-j} + \sum_{k=1}^K \delta_{k} Z_{k} + \epsilon_{t} + \left\lbrace
#'     \begin{array}{ll}
#' \phi_{1,0} + \sum_{i \in I_1} \phi_{1,i} X_{t-i}   & \mathrm{if } X_{t-d} \leq \mathrm{thd} \\\\\\
#'  &\\\\\\
#' \phi_{2,0} + \sum_{i \in I_2} \phi_{2,i} X_{t-i}   & \mathrm{if } X_{t-d} > \mathrm{thd}
#' \end{array}
#' \right. }{X[t] = \phi[0] + \Sigma_{h in I}  \phi[h] X[t-h] + \Sigma_{j = 1,..,q} \theta[j] \epsilon[t-j] + \Sigma_{j = 1,..,Q} \Theta[j] \epsilon[t-js] + \Sigma_{k = 1,..,K} \delta[k] Z[k] + \epsilon[t] +
#' + \phi[1,0] + \Sigma_{i  in I_1} \phi[1,i] X[t-i] --  if X[t-d] <= thd
#' + \phi[2,0] + \Sigma_{i in I_2} \phi[2,i] X[t-i]  --  if X[t-d] > thd}
#'
#' where \mjeqn{\phi_h}{\phi[h]} are the common AR parameters and \mjseqn{h} ranges in \code{I = ar.lags}. \mjeqn{\theta_j}{\theta[j]} are the common MA parameters and \mjeqn{j = 1,\dots,q}{j = 1,...,q}
#' (\code{q = ma.ord}), \mjeqn{\Theta_l}{\Theta[l]} are the common seasonal MA parameters and \mjeqn{l = 1,\dots,Q}{l = 1,...,Q} (\code{Q = sma.ord})
#' \mjeqn{\delta_k}{\delta[k]} are the parameters for the covariates. Finally, \mjeqn{\phi_{1,i}}{\phi[1,i]} and \mjeqn{\phi_{2,i}}{\phi[2,i]} are the TAR parameters
#' for the lower and upper regime, respectively and \code{I1 = tar1.lags} \code{I2 = tar2.lags} are the vector of TAR lags.
#'
#' @return
#'   A list of class \code{TARMA} with components:
#' \itemize{
#'  \item \code{fit}    - The output of the fit. It is a \code{arima} object.
#'  \item \code{aic}    - Value of the AIC for the minimised least squares criterion over the threshold range.
#'  \item \code{bic}    - Value of the BIC for the minimised least squares criterion over the threshold range.
#'  \item \code{aic.v}  - Vector of values of the AIC over the threshold range.
#'  \item \code{thd.range} - Vector of values of the threshold range.
#'  \item \code{d}      - Delay parameter.
#'  \item \code{thd}    - Estimated threshold.
#'  \item \code{phi1}   - Estimated AR parameters for the lower regime.
#'  \item \code{phi2}   - Estimated AR parameters for the upper regime.
#'  \item \code{theta1} - Estimated MA parameters for the lower regime.
#'  \item \code{theta2} - Estimated MA parameters for the upper regime.
#'  \item \code{delta}  - Estimated parameters for the covariates \code{x.reg}.
#'  \item \code{tlag1}  - TAR lags for the lower regime
#'  \item \code{tlag2}  - TAR lags for the upper regime
#'  \item \code{mlag1}  - TMA lags for the lower regime
#'  \item \code{mlag2}  - TMA lags for the upper regime
#'  \item \code{arlag}  -  Same as the input slot \code{ar.lags}
#'  \item \code{include.int} -  Same as the input slot \code{include.int}
#'  \item \code{se}     - Standard errors for the parameters. Note that they are computed conditionally upon the threshold so that they are generally smaller than the true ones.
#'  \item \code{rss}    - Minimised residual sum of squares.
#'  \item \code{method} - Estimation method.
#'  \item \code{call}   - The matched call.
#'}
#' @author Simone Giannerini, \email{simone.giannerini@@unibo.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Gia21}{tseriesTARMA}
#' * \insertRef{Cha19}{tseriesTARMA}
#'
#' @importFrom stats is.ts median quantile arima lag ts.intersect
#' @export
#' @seealso \code{\link{TARMA.fit}} for Least Square estimation of full subset TARMA
#' models. \code{\link{print.TARMA}} for print methods for \code{TARMA} fits.
#' \code{\link{predict.TARMA}} for prediction and forecasting.
#' @examples
#' ## a TARMA(1,1,1,1)
#' set.seed(127)
#' x    <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0,0.8), theta1=0.5, theta2=0.5, d=1, thd=0.2)
#' fit1 <- TARMA.fit2(x, tar1.lags=1, tar2.lags=1, ma.ord=1, d=1)
#'
#' \donttest{
#' ## Showcase of the fit with covariates ---
#' ## simulates from a TARMA(3,3,1,1) model with common MA parameter
#' ## and common AR(1) and AR(2) parameters. Only the lag 3 parameter varies across regimes
#' set.seed(212)
#' n <- 300
#' x <- TARMA.sim(n=n, phi1=c(0.5,0.3,0.2,0.4), phi2=c(0.5,0.3,0.2,-0.2), theta1=0.4, theta2=0.4,
#'      d=1, thd=0.2, s1=1, s2=1)
#'
#' ## FIT 1: estimates lags 1,2,3 as threshold lags ---
#' fit1 <- TARMA.fit2(x, ma.ord=1, tar1.lags=c(1,2,3), tar2.lags=c(1,2,3), d=1)
#'
#' ## FIT 2: estimates lags 1 and 2 as fixed AR and lag 3 as the threshold lag
#' fit2 <- TARMA.fit2(x, ma.ord=1, tar1.lags=c(3),  tar2.lags=c(3), ar.lags=c(1,2), d=1)
#'
#' ## FIT 3: creates lag 1 and 2 and fits them as covariates ---
#' z1   <- lag(x,-1)
#' z2   <- lag(x,-2)
#' fit3 <- TARMA.fit2(x, ma.ord=1,  tar1.lags=c(3), tar2.lags=c(3), x.reg=ts.intersect(z1,z2), d=1)
#'
#' ## FIT 4: estimates lag 1 as a covariate, lag 2 as fixed AR and lag 3 as the threshold lag
#' fit4 <- TARMA.fit2(x, ma.ord = 1,  tar1.lags=c(3), tar2.lags=c(3), x.reg=z1, ar.lags=2, d=1)
#'
#' }
## '***************************************************************************

TARMA.fit2 <- function(x, ar.lags = NULL, tar1.lags = c(1), tar2.lags = c(1), ma.ord = 1, sma.ord = 0L,
period = NA, estimate.thd = TRUE, threshold, d = 1, pa = 0.25, pb = 0.75,
thd.var = NULL, include.int = TRUE, x.reg = NULL, optim.control = list(), ...){

    if(length(intersect(ar.lags,tar1.lags))>0) stop('AR lags and TAR lags must not have common elements')
    if(length(intersect(ar.lags,tar2.lags))>0) stop('AR lags and TAR lags must not have common elements')
    if(!is.ts(x)) x <- ts(x)
    xtsp    <- tsp(x)
    n       <- length(x)
    n.ar    <- length(ar.lags)
    n.tar1  <- length(tar1.lags)
    n.tar2  <- length(tar2.lags)
    t.lags  <- c(ar.lags,tar1.lags,tar2.lags)
    if(any(t.lags<=0)) stop('lags must be positive integers')
    n.tot   <- length(t.lags)
    if(is.null(thd.var)){
        thd.var <- lag(x,-d)
    }else{
        d <- NA
        if(!is.ts(thd.var)) stop('thd.var needs to be a time series. Check the time alignment with x')
    }
    xx    <- ts.intersect(x,thd.var)
    if(n.ar>0){
         arnames  <- paste('ar',ar.lags,sep='')
    } else {
         arnames <- NULL
    }
    for(k in 1:n.tot){
        xx <- ts.intersect(xx,lag(x,-t.lags[k]))
    }
    # add regressors, if there are any
    if(!is.null(x.reg)){
      if(NROW(x.reg)!=n) warning('The regressors and the input time series x have different lengths')  
      if(!is.ts(x.reg)){
        x.reg <- ts(x.reg,start=xtsp[1],frequency = xtsp[3])
      }
        n.reg <- NCOL(x.reg)
        xx    <- ts.intersect(xx,x.reg)
    }else{
      n.reg <-0
    }
    xth   <- xx[,2]          # threshold variable
    neff  <- length(xth)     # effective sample size
    xreg  <- xx[,-(1:2),drop=FALSE]   # ar+tar+regression lags
    nreg  <- NCOL(xreg)      # total number of ar+tar+regression parameters
    if(n.ar>0){
        xar <- xreg[,1:n.ar,drop=FALSE]
    } else {
        xar <- NULL
    }
    xtar1 <- xreg[,(n.ar+1):(n.ar+n.tar1)] # TAR1 regressors
    xtar2 <- xreg[,(n.ar+n.tar1+1):(n.ar+n.tar1+n.tar2)] # TAR2 regressors
    if(!is.null(x.reg)){
        x.reg  <- xreg[,(n.tar1+n.tar2+n.ar+1):nreg]    # covariates
        cn.reg <- paste('Z',1:n.reg)
    }
    if(include.int){ # adds the intercept to the tar model
        xtar1    <- cbind(rep(1,neff),xtar1)
        xtar2    <- cbind(rep(1,neff),xtar2)
        ntar1   <- c('int.1',paste('ar1',tar1.lags,sep='.'))
        ntar2   <- c('int.2',paste('ar2',tar2.lags,sep='.'))
    }else{    # adds an overall intercept
        xar     <- cbind(rep(1,neff),xar)
        arnames <- c('int',arnames)
        ntar1   <- paste('ar1',tar1.lags,sep='.')
        ntar2   <- paste('ar2',tar2.lags,sep='.')
    }
    tarnames <- c(ntar1,ntar2)
    if(!is.null(xar)){
        xar            <- as.matrix(xar)
        colnames(xar)  <- arnames
    }
    xtar1  <- as.matrix(xtar1)
    xtar2  <- as.matrix(xtar2)
#    colnames(xtar) <- tarnames
    n.tarp <- length(tarnames)  # number of TAR parameters
    xnames <- c(arnames,tarnames)

    if(estimate.thd){
        a         <- ceiling((neff-1)*pa)
        b         <- floor((neff-1)*pb)
        thd.v     <- sort(xth)
        thd.range <- thd.v[a:b]
        nr        <- length(thd.range)
        rss.g     <- Inf
        rssv      <- rep(NA,nr)
        ind <- trunc(quantile(1:nr,prob=seq(0,1,by=0.1)))
        for(i in 1:nr){
            thd     <- thd.range[i]
            Il      <- (xth<= thd)
            Iu      <- (xth>  thd)
            xtar.h  <- matrix(c(xtar1*Il,xtar2*Iu),nrow=neff,ncol=n.tarp)
            xreg.h  <- cbind(xar,xtar.h,x.reg)
            fit.tar <- tryCatch(stats::arima(xx[,1],order=c(0L,0L,ma.ord),
                transform.pars = FALSE, include.mean=FALSE, xreg=xreg.h,
                seasonal = list(order = c(0L, 0L, sma.ord), period = period),...)
                , error=function(x){NULL})
            if(!is.null(fit.tar)){
                resi    <- fit.tar$residuals  #\hat{e_t}
                rss.v   <- fit.tar$aic
                rssv[i] <- rss.v
            }else{
                rss.v <- Inf
            }
            if(rss.v < rss.g){
                rss.g  <- rss.v
                resi.g <- resi
                fit.g  <- fit.tar
                igood  <- i
            }
        }
        thd       <- thd.range[igood]
    }else{
        thd       <- threshold
        rssv      <- NA
        thd.range <- NA
    }
    Il        <- (xth<= thd)
    Iu        <- (xth>  thd)
    xtar.h  <- matrix(c(xtar1*Il,xtar2*Iu),nrow=neff,ncol=n.tarp)
    xreg.h  <- cbind(xar,xtar.h,x.reg)
    if(is.null(x.reg)){
        colnames(xreg.h) <- xnames
    }else{
        colnames(xreg.h) <- c(xnames,cn.reg)
    }
    fit.g <- stats::arima(xx[,1],order=c(0L,0L,ma.ord),
    transform.pars = FALSE, include.mean=FALSE, xreg=xreg.h,
    seasonal = list(order = c(0L, 0L, sma.ord), period = period),...)
    pars   <- fit.g$coef
    if(is.na(period)){
        period   <- 0
        sma.lags <- 0
    }else{
        sma.lags <- (1:sma.ord)*period
    }
    nma    <- max(ma.ord,sma.lags)
    theta1 <- rep(0L,nma)
    names(theta1) <- 1:nma
    theta2 <- theta1
    if(nma>0){
        theta1[1:ma.ord] <- pars[1:ma.ord]
        theta2[1:ma.ord] <- pars[1:ma.ord]
        pars   <- pars[-(1:ma.ord)]
        if(sma.ord>0){
            for(i in 1:sma.ord){
                theta1[sma.lags[i]] <- pars[i]
                theta2[sma.lags[i]] <- pars[i]
            }
            pars <- pars[-(1:sma.ord)]
        }
    }
    mlag     <- which(theta1!=0)
    if(!include.int){ # common
        int.1 <- int.2 <- pars[1]
        pars  <- pars[-1]
    }else{
        int.1 <- pars['int.1']
        int.2 <- pars['int.2']
        pars  <- pars[-c(which(names(pars)=='int.1'),which(names(pars)=='int.2'))]
    }
    if(n.ar>0){
        p.ar   <- pars[1:n.ar]
        pars   <- pars[-(1:n.ar)]
        ind.ar <- ar.lags+1
    }
    p.1      <- pars[1:n.tar1]
    p.2      <- pars[(n.tar1+1):(n.tar1+n.tar2)]
    pars     <- pars[-(1:(n.tar1+n.tar2))]
    ind.tar1 <- tar1.lags+1
    ind.tar2 <- tar2.lags+1
    phi1     <- rep(0,max(t.lags)+1)
    names(phi1) <- 0:max(t.lags)
    phi2        <- phi1
    phi1[1]     <- int.1
    phi2[1]     <- int.2
    if(n.ar>0){
        phi1[ind.ar]  <- p.ar
        phi2[ind.ar]  <- p.ar
    }
    phi1[ind.tar1] <- p.1
    phi2[ind.tar2] <- p.2
    p.Z     <- pars
    # computes AIC and BIC
    Ir      <- Il
    n1      <- sum(Ir)
    n2      <- sum(1-Ir)
    np1     <- n.tar1+1
    np2     <- n.tar2+1
    npars   <- n.ar+length(mlag)+np1+np2+1-(1-include.int) # subtracts 1 in case of common intercept
    eps     <- residuals(fit.g)
    s2hat.1 <- sum(Ir*eps^2)/n1
    s2hat.2 <- sum((1-Ir)*eps^2)/n2
    aic     <- n1*log(s2hat.1) + n2*log(s2hat.2) + neff*(1+log(2*pi)) + 2*npars
    bic     <- n1*log(s2hat.1) + n2*log(s2hat.2) + neff*(1+log(2*pi)) + log(neff)*npars

    res  <- list(fit=fit.g, aic=aic, aic.v=rssv, bic=bic, thd.range=thd.range, thd=thd,
        d = d, phi1=phi1, phi2=phi2, theta1=theta1, theta2=theta2, delta=p.Z,
        tlag1=tar1.lags, tlag2=tar2.lags,mlag1=mlag,mlag2=mlag,arlag=ar.lags,include.int=include.int,
        se=sqrt(diag(fit.g$var.coef)), rss=sum(fit.g$residuals^2,na.rm=TRUE),
        method='MLE', call = match.call())
    class(res) <- 'TARMA'
    return(res)
}
## ***************************************************************************
