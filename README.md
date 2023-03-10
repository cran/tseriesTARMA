
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tseriesTARMA

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/tseriesTARMA)](https://CRAN.R-project.org/package=tseriesTARMA)
[![CRAN_download](http://cranlogs.r-pkg.org/badges/tseriesTARMA)](https://cran.r-project.org/package=tseriesTARMA)
[![CRAN_download_total](http://cranlogs.r-pkg.org/badges/grand-total/tseriesTARMA)](https://cran.r-project.org/package=tseriesTARMA)
<!-- badges: end -->

# Analysis of Nonlinear Time Series through Threshold Autoregressive Moving Average Models (TARMA) models

It provides advanced functions for:

- TARMA model fitting and forecasting:
  - Least Squares fitting of a full subset TARMA model, including robust
    LS fitting.
  - Maximum Likelihood fitting of a subset TARMA model with common MA
    parts and possible covariates.
- TARMA testing for threshold type nonlinearity:
  - Tests for AR vs TAR (asymptotic, bootstrap, wild bootstrap)
  - Tests for ARMA vs TARMA with both i.i.d. errors and GARCH errors.
- Unit-root testing against a stationary TARMA model

## Installation

``` r
install.packages("tseriesTARMA")
```

## Authors

- [Simone Giannerini, University of
  Bologna](https://www.unibo.it/sitoweb/simone.giannerini/en)
- [Greta Goracci, Free University of
  Bolzano/Bozen](https://www.unibz.it/it/faculties/economics-management/academic-staff/person/46136-greta-goracci)

## References

- Giannerini, Goracci, and Rahbek (2023)
- Giannerini, Goracci, and Rahbek (2022)
- Giannerini and Goracci (2021)
- Goracci et al. (2021)
- Goracci et al. (2023)
- K.-S. Chan and Goracci (2019)
- K.-S. Chan et al. (2024)

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Cha19" class="csl-entry">

Chan, K. -S., and G. Goracci. 2019. βOn the Ergodicity of First-Order
Threshold Autoregressive Moving-Average Processes.β *J. Time Series
Anal.* 40 (2): 256β64.

</div>

<div id="ref-Cha24" class="csl-entry">

Chan, K.-S., S. Giannerini, G. Goracci, and H. Tong. 2024. βTesting for
Threshold Regulation in Presence of Measurement Error.β *Statistica
Sinica* 34 (3). <https://doi.org/10.5705/ss.202022.0125>.

</div>

<div id="ref-Gia21" class="csl-entry">

Giannerini, S., and G. Goracci. 2021. βEstimating and Forecasting with
TARMA Models.β University of Bologna.

</div>

<div id="ref-Gia22" class="csl-entry">

Giannerini, S., G. Goracci, and A. Rahbek. 2022. βThe Validity of
Bootstrap Testing in the Threshold Framework.β arXiv.
<https://doi.org/10.48550/ARXIV.2201.00028>.

</div>

<div id="ref-Gia23" class="csl-entry">

βββ. 2023. βThe Validity of Bootstrap Testing in the Threshold
Framework.β *Journal of Econometrics*, forthcoming.
<https://doi.org/10.1016/j.jeconom.2023.01.004>.

</div>

<div id="ref-Gor21" class="csl-entry">

Goracci, G., S. Giannerini, K.-S. Chan, and H. Tong. 2021. βTesting for
Threshold Effects in the TARMA Framework.β University of Bologna, Free
University of Bolzano, University of Iowa, London School of Economics.
<https://arxiv.org/abs/2103.13977>.

</div>

<div id="ref-Gor23" class="csl-entry">

βββ. 2023. βTesting for Threshold Effects in the TARMA Framework.β
*Statistica Sinica* 33 (3): to appear.
<https://doi.org/10.5705/ss.202021.0120>.

</div>

</div>
