#' @title Methods for TARMA fits
#'
#' @param x       A \code{TARMA} fit.
#' @param object  A \code{TARMA} fit.
#' @param digits Number of decimal digits for the output.
#' @param se Logical. if \code{TRUE} (the default) prints also the standard errors.
#' @param \dots Further parameters.
#' @return
#'   No return value, called for side effects
#' @method print TARMA
#' @export
#' @seealso \code{\link{TARMA.fit}} and \code{\link{TARMA.fit2}} for TARMA modelling, \code{\link{plot.tsfit}} for plotting TARMA fits and forecasts.


print.TARMA <-
    function (x, digits = max(3L, getOption("digits") - 3L), se = TRUE, ...){
    if(x$method=='MLE'){
        mask      <- x$fit$mask
        i.int     <- x$include.int
        n.ma      <- length(x$mlag1)       # number of common MA+SMA parameters
        n.ar      <- length(x$arlag)       # number of common AR parameters
        n.tar1    <- length(x$tlag1)+i.int # number of TAR parameters (lower)
        n.tar2    <- length(x$tlag2)+i.int # number of TAR parameters (upper)
        nptot     <- length(x$se)          # total number of parameters
        coefm     <- round(x$fit$coef, digits = digits)    
        coef.arma <- coefm[mask][1:(n.ma+n.ar+(1-i.int))]
        coef.l    <- coefm[mask][((n.ma+n.ar+(1-i.int))+1):((n.ma+n.ar+(1-i.int))+n.tar1)]
        coef.u    <- coefm[mask][((n.ma+n.ar+(1-i.int))+n.tar1+1):((n.ma+n.ar+(1-i.int))+n.tar1+n.tar2)]
        nxreg     <- length(x$delta)
        if(nxreg){
            coef.reg <- x$delta
        }
        if(se){
            ses      <- rep_len(0, nptot)
            ses[x$fit$mask] <- round(x$se, digits = digits)
            ses.arma <- ses[mask][1:(n.ma+n.ar+(1-i.int))]
            ses.l    <- ses[mask][(n.ma+n.ar+(1-i.int)+1):(n.ma+n.ar+(1-i.int)+n.tar1)]
            ses.u    <- ses[mask][(n.ma+n.ar+(1-i.int)+n.tar1+1):(n.ma+n.ar+(1-i.int)+n.tar1+n.tar2)]
            if(nxreg){
                ses.reg <- ses[mask][(n.ma+n.ar+(1-i.int)+n.tar1+n.tar2+1):nptot]
            }
            coef.arma <- matrix(coef.arma, 1L, dimnames = list(NULL, names(coef.arma)))
            coef.arma <- rbind(coef.arma, s.e. = ses.arma)
            coef.l    <- matrix(coef.l, 1L, dimnames = list(NULL, names(coef.l)))
            coef.l    <- rbind(coef.l, s.e. = ses.l)
            coef.u    <- matrix(coef.u, 1L, dimnames = list(NULL, names(coef.u)))
            coef.u    <- rbind(coef.u, s.e. = ses.u)
            if(nxreg){
                coef.reg    <- matrix(coef.reg, 1L, dimnames = list(NULL, names(coef.reg)))
                coef.reg    <- rbind(coef.reg, s.e. = ses.reg)
            }
        }
        cat('-----------------------------------------------------------------\n')
        cat("--- Maximum Likelihood fit of the TARMA model:\n")
        cat('-----------------------------------------------------------------\n')
        cat("\nCall:", deparse(x$call, width.cutoff = 75L), "", sep = "\n")
        cat('-----------------------------------------------------------------\n')
        cat(" - Common ARMA coefficients:\n\n")
        print.default(coef.arma, print.gap = 2)
        if(nxreg){
        cat('-----------------------------------------------------------------\n')
        cat(" - Covariates:\n")
        print.default(coef.reg, print.gap = 2, digits = digits)
        }
        cat('-----------------------------------------------------------------\n')
        cat(" - Lower regime: ")
        cat("TAR lags: ", x$tlag1, "\n\n", sep = " ")
        cat("Coefficients:\n")
        print.default(coef.l, print.gap = 2)
        cat('-----------------------------------------------------------------\n')
        cat(" - Upper regime: ")
        cat("TAR lags: ", x$tlag2, "\n\n", sep = " ")
        cat("Coefficients:\n")
        print.default(coef.u, print.gap = 2)
    }else{
        np1 <- length(x$tlag1)+1
        np2 <- length(x$tlag2)+1
        nq1 <- length(x$mlag1)
        nq2 <- length(x$mlag2)
        dum <- ifelse(x$method=='robust','Robust M','Least Squares')
        cat('-----------------------------------------------------------------\n')
        cat('--- ',dum,'fit of the TARMA model:\n')
        cat('-----------------------------------------------------------------\n')
        cat('\nCall:', deparse(x$call, width.cutoff = 75L), '', sep = "\n")
        cat('Method: ', x$method, '  - Innovation density: ', x$innov, '\n', sep = '')
        if (length(x$fit$coef)){
            coef.l <- round(c(x$fit$coef[1:np1],x$fit$coef[(np1+np2+1):(np1+np2+nq1)]), digits = digits)
            coef.u <- round(c(x$fit$coef[(np1+1):(np1+np2)],x$fit$coef[(np1+np2+nq1+1):(np1+np2+nq1+nq2)]), digits = digits)
            if (se) {
                ses.l  <- rep.int(0, length(coef.l))
                ses.l  <- round(c(x$se[1:np1],x$se[(np1+np2+1):(np1+np2+nq1)]), digits = digits)
                coef.l <- matrix(coef.l, 1L, dimnames = list(NULL, names(coef.l)))
                coef.l <- rbind(coef.l, s.e. = ses.l)
                ses.u  <- rep.int(0, length(coef.u))
                ses.u  <- round(c(x$se[(np1+1):(np1+np2)],x$se[(np1+np2+nq1+1):(np1+np2+nq1+nq2)]), digits = digits)
                coef.u <- matrix(coef.u, 1L, dimnames = list(NULL, names(coef.u)))
                coef.u <- rbind(coef.u, s.e. = ses.u)
            }
        cat('-----------------------------------------------------------------\n')
            cat(" - Lower regime:\n\n")
            cat("TAR lags: ", x$tlag1, "\n", sep = " ")
            cat("TMA lags: ", x$mlag1, "\n", sep = " ")
            cat("\nCoefficients:\n")
            print.default(coef.l, print.gap = 2)
        cat('-----------------------------------------------------------------\n')
            cat(" - Upper regime:\n\n")
            cat("TAR lags: ", x$tlag2, "\n", sep = " ")
            cat("TMA lags: ", x$mlag2, "\n", sep = " ")
            cat("\nCoefficients:\n")
            print.default(coef.u, print.gap = 2)
        }
    }
    aicdum   <- NULL
    alphadum <- NULL
    if(!is.null(x$aic)){aicdum <- paste("AIC = ",format(round(x$aic, 2L)),sep='')}
    if(!is.null(x$alpha)&&x$alpha>0){alphadum <- paste("alpha = ",format(round(x$alpha, 2L)), "; qu = ",toString(format(round(x$qu, 2L))),sep='')}
        cat('-----------------------------------------------------------------\n')
    cat("\nd = ",x$d,";  threshold = ", format(round(x$thd, 2L)),
        "; ", alphadum, "\n", sep = "")
    cat('-----------------------------------------------------------------\n\n')
    cat("sigma^2 estimated as ", format(x$fit$sigma2, digits = digits),
        "; RSS = ", format(round(x$rss, 2L)),"; ", aicdum, "\n", sep = "")
    cat('-----------------------------------------------------------------\n')
    if(!is.null(x$innovpar)){
      cat("Innovations parameter estimated as", format(x$innovpar, digits = digits),'\n\n')
    }
    invisible(x)
}

#' @rdname print.TARMA
#' @method coef TARMA
#' @export
#'
coef.TARMA <- function (object, ...) object$fit$coef

#' @rdname print.TARMA
#' @method vcov TARMA
#' @export
#'
vcov.TARMA <- function (object, ...) object$fit$var.coef

#' @rdname print.TARMA
#' @method residuals TARMA
#' @export
#'
residuals.TARMA <- function (object, ...) object$fit$residuals
