#' Simulation of a two-regime \code{TARMA(p1,p2,q1,q2)} process
#'
#' @description \loadmathjax
#' Simulates from the following two-regime \code{TARMA(p1,p2,q1,q2)} process: \cr
#'
#' \mjdeqn{X_{t} = \left\lbrace
#' \begin{array}{ll}
#'  \phi_{1,0} + \sum_{i=1}^{p_1} \phi_{1,i} X_{t-i} + \sum_{j=1}^{q_1} \theta_{1,j} \varepsilon_{t-j} + \varepsilon_{t}, & \quad\mathrm{if}\quad X_{t-d} \leq \mathrm{thd} \\\\\\
#'  &\\\\\\
#'  \phi_{2,0} + \sum_{i=1}^{p_2} \phi_{2,i} X_{t-i} + \sum_{j=1}^{q_2} \theta_{2,j} \varepsilon_{t-j} + \varepsilon_{t}, & \quad\mathrm{if}\quad X_{t-d} > \mathrm{thd}
#' \end{array}
#' \right.
#' }{X[t] =
#'   = \phi[1,0] + \phi[1,1]X[t-1] + ... + \phi[1,p1]X[t-p1] + \theta[1,1]\epsilon[t-1] + ... +  \theta[1,q]\epsilon[t-q1] + \epsilon[t] --   if  X[t-d] <=  thd
#'
#'   = \phi[2,0] + \phi[2,1]X[t-1] + ... + \phi[2,p2]X[t-p2] + \theta[2,1]\epsilon[t-1] + ... +  \theta[2,q]\epsilon[t-q2]  + \epsilon[t] --  if  X[t-d] >  thd}
#'
#'
#' @param n      Length of the series.
#' @param phi1   Vector of \code{p1+1} Autoregressive parameters of the lower regime.
#'               The first element is the intercept.
#' @param phi2   Vector of \code{p2+1} Autoregressive parameters of the upper regime.
#'               The first element is the intercept.
#' @param theta1 Vector of \code{q1} Moving Average parameters of the lower regime.
#' @param theta2 Vector of \code{q2} Moving Average parameters of the upper regime.
#' @param d      Delay parameter. Defaults to \code{1}.
#' @param thd    Threshold parameter. Defaults to \code{0}.
#' @param s1     Innovation variance for the lower regime. Defaults to \code{1}.
#' @param s2     Innovation variance for the upper regime. Defaults to \code{1}.
#' @param rand.gen Optional: a function to generate the innovations. Defaults to \code{rnorm}.
#' @param innov   Optional: a time series of innovations. If not provided, \code{rand.gen} is used.
#' @param n.start Length of the burn-in period. Defaults to \code{500}.
#' @param start.innov   Optional: a time series of innovations for the burn-in period.
#' @param xstart  Initial condition as a named list: \cr
#'                \code{$ar}: AR part length \code{k = max(p1,p2,d),  X[k], X[k-1], ... ,X[1]}; \cr
#'                \code{$ma}: MA part length \code{q = ma.ord,   e[q], ... , e[1]}.
#' @param \dots  Additional arguments for \code{rand.gen}.
#'
#' @details Note that the parameters are not checked for ergodicity.
#'
#' @return A time series object of class \code{ts} generated from the above model.
#'
#' @examples
#' ## a TARMA(1,1,1,1) model
#' set.seed(123)
#' x <- TARMA.sim(n=100, phi1=c(0.5,-0.5), phi2=c(0.0,0.8), theta1=-0.5, theta2=0.5, d=1, thd=0.2)
#'
#' ## a TARMA(1,2,1,1) model
#' x <- TARMA.sim(n=100,phi1=c(0.5,-0.5,0),phi2=c(0,0.5,0.3),theta1=-0.5,theta2=0.5,d=1,thd=0.2)
#'
#' @author Simone Giannerini, \email{simone.giannerini@@unibo.it}
#' @author Greta Goracci, \email{greta.goracci@@unibz.it}
#' @references
#' * \insertRef{Gia21}{tseriesTARMA}
#'
#' @importFrom stats ts rnorm
#' @export
#'

## ****************************************************************************

TARMA.sim  <- function(n, phi1, phi2, theta1, theta2, d=1, thd=0, s1=1, s2=1, rand.gen = rnorm,
innov = rand.gen(n, ...), n.start = 500, xstart, start.innov = rand.gen(n.start, ...),...){

    if(n.start<0) stop('n.start must be greater than 0')
#    if((length(theta1)!=length(theta2)))
#    stop('the parameters vectors theta must have the same length')

    tar1.ord <- length(phi1)-1
    tar2.ord <- length(phi2)-1
    tma1.ord <- length(theta1)
    tma2.ord <- length(theta2)
    k        <- max(tar1.ord,tar2.ord,d,tma1.ord,tma2.ord)
    if(missing(xstart)) xstart <- list(ar=rep(0,k),ma=rep(0,k))
    mlag   <- max(k,d)
    n.d    <- mlag + n.start # number of discarded observations
    if(n.start==0){
        start.innov <- rep(0,mlag)
    } else{
        start.innov <- c(rand.gen(mlag),start.innov)
    }
    e      <- c(start.innov, innov[1L:n])
    ntot   <- length(e)
    x      <- double(ntot)
    x[mlag:(mlag-k+1)] <- xstart$ar;
    e[mlag:(mlag-k+1)] <- xstart$ma
    for(i in (mlag+1):ntot){
        if(x[i-d] <= thd){
            x[i] <- (phi1%*%c(1,x[(i-1):(i-tar1.ord)]))[1] + s1*e[i] + theta1%*%c(s1*e[(i-1):(i-tma1.ord)])
        }else{
            x[i] <- (phi2%*%c(1,x[(i-1):(i-tar2.ord)]))[1] + s2*e[i] + theta2%*%c(s2*e[(i-1):(i-tma2.ord)])
        }
    }
    x <- x[-(1:n.d)]
    return(ts(x));
}
## ***************************************************************************
