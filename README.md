[![CRAN_version](http://www.r-pkg.org/badges/version/tseriesTARMA)](https://cran.r-project.org/package=tseriesTARMA)
[![CRAN_download](http://cranlogs.r-pkg.org/badges/tseriesTARMA)](https://cran.r-project.org/package=tseriesTARMA)
[![CRAN_download_total](http://cranlogs.r-pkg.org/badges/grand-total/tseriesTARMA)](https://cran.r-project.org/package=tseriesTARMA)

***
# tseriesTARMA

# Analysis of Nonlinear Time Series through Threshold Autoregressive Moving Average Models (TARMA) models 

It provides advanced functions for:

 * TARMA model fitting and forecasting:
    - Least Squares fitting of a full subset TARMA model, including robust LS fitting.
    - Maximum Likelihood fitting of a subset TARMA model with common MA parts and possible covariates.
    
 * TARMA testing for threshold type nonlinearity:
    - Tests for AR vs TAR (asymptotic, bootstrap, wild bootstrap)
    - Tests for ARMA vs TARMA with both i.i.d. errors and GARCH errors.
    
 * Unit-root testing against a stationary TARMA model


## Install
```r
install.packages("tseriesTARMA")
```

## Authors
 - [Simone Giannerini, University of Bologna](https://www.unibo.it/sitoweb/simone.giannerini/en) 
 - [Greta Goracci, Free University of Bolzano/Bozen](https://www.unibz.it/it/faculties/economics-management/academic-staff/person/46136-greta-goracci)

## References
 * \insertRef{Gia23}{tseriesTARMA}
 * \insertRef{Gia22}{tseriesTARMA}
 * \insertRef{Gia21}{tseriesTARMA}
 * \insertRef{Gor21}{tseriesTARMA}
 * \insertRef{Gor23}{tseriesTARMA}
 * \insertRef{Cha19}{tseriesTARMA}
 * \insertRef{Cha24}{tseriesTARMA}