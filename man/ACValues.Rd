% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ACValues.R
\docType{data}
\name{ACValues}
\alias{ACValues}
\title{Andrews Tabulated Critical Values}
\format{
\subsection{\code{ACValues}}{

A matrix with 13 rows and 62 columns. The first two columns contain the parameters, whereas
the remaining 60 columns contain 3 critical values (at level 10\%, 5\%, and 1\%) for each value of the
parameter \mjeqn{p}{p}:
\describe{
\item{\code{pi}}{\mjeqn{\pi_0}{\pi[0]} in the paper, corresponds to \code{pa}.}
\item{\code{lambda}}{\mjeqn{\lambda}{\lambda} in the paper, not relevant here.}
\item{\code{1-90}}{\mjeqn{p=1}{p=1}, critical level 10\%.}
\item{\code{1-95}}{\mjeqn{p=1}{p=1}, critical level 5\%.}
\item{\code{1-99}}{\mjeqn{p=1}{p=1}, critical level 1\%.}
...
\item{\code{20-90}}{\mjeqn{p=20}{p=20}, critical level 10\%.}
\item{\code{20-95}}{\mjeqn{p=20}{p=20}, critical level 5\%.}
\item{\code{20-99}}{\mjeqn{p=20}{p=20}, critical level 1\%.}
}
Note that \mjeqn{p}{p} are the degrees of freedom and correspond to the total number of tested parameter
in the above tests.
}
}
\source{
\insertRef{And03}{tseriesTARMA}
}
\usage{
ACValues
}
\description{
\loadmathjax
The data is taken from Table 1 of \insertCite{And03}{tseriesTARMA}, which provides asymptotic critical
values for sup Wald, LM, and LR tests for parameter instability. Critical values are given for
degrees of freedom \mjeqn{p=1,\dots,20}{p = 1, ..., 20}.
They can be used as asymptotic critical values for the following tests:
\itemize{
\item \code{\link{TAR.test}}
\item \code{\link{TARMA.test}}
\item \code{\link{TARMAGARCH.test}}
}

Provided \code{pb = 1- pa}.
}
\keyword{datasets}
