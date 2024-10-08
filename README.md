
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
    fitting based on M-estimators.
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
  Udine](https://www.simonegiannerini.net)
- [Greta Goracci, Free University of
  Bolzano/Bozen](https://www.unibz.it/it/faculties/economics-management/academic-staff/person/46136-greta-goracci)

## References

- Goracci et al. (2025)
- Angelini et al. (2023)
- Giannerini, Goracci, and Rahbek (2024)
- Goracci, Ferrari, et al. (2023)
- Giannerini, Goracci, and Rahbek (2022)
- Giannerini and Goracci (2021)
- Goracci et al. (2021)
- Goracci, Giannerini, et al. (2023)
- K.-S. Chan and Goracci (2019)
- K.-S. Chan et al. (2024)

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-Ang23" class="csl-entry">

Angelini, F., M. Castellani, S. Giannerini, and G. Goracci. 2023.
“Testing for Threshold Effects in Presence of Heteroskedasticity and
Measurement Error with an Application to Italian Strikes.” University of
Bologna; Free University of Bolzano. <https://arxiv.org/abs/2308.00444>.

</div>

<div id="ref-Cha19" class="csl-entry">

Chan, K. -S., and G. Goracci. 2019. “On the Ergodicity of First-Order
Threshold Autoregressive Moving-Average Processes.” *J. Time Series
Anal.* 40 (2): 256–64.

</div>

<div id="ref-Cha24" class="csl-entry">

Chan, K.-S., S. Giannerini, G. Goracci, and H. Tong. 2024. “Testing for
Threshold Regulation in Presence of Measurement Error.” *Statistica
Sinica* 34 (3): 1413–34. <https://doi.org/10.5705/ss.202022.0125>.

</div>

<div id="ref-Gia21" class="csl-entry">

Giannerini, S., and G. Goracci. 2021. “Estimating and Forecasting with
TARMA Models.” University of Bologna.

</div>

<div id="ref-Gia22" class="csl-entry">

Giannerini, S., G. Goracci, and A. Rahbek. 2022. “The Validity of
Bootstrap Testing in the Threshold Framework.” arXiv.
<https://doi.org/10.48550/ARXIV.2201.00028>.

</div>

<div id="ref-Gia23" class="csl-entry">

———. 2024. “The Validity of Bootstrap Testing in the Threshold
Framework.” *Journal of Econometrics* 239 (1): 105379.
<https://doi.org/10.1016/j.jeconom.2023.01.004>.

</div>

<div id="ref-Gor23b" class="csl-entry">

Goracci, G., D. Ferrari, S. Giannerini, and F. Ravazzolo. 2023. “Robust
Estimation for Threshold Autoregressive Moving-Average Models.” Free
University of Bolzano, University of Bologna.
<https://doi.org/10.48550/ARXIV.2211.08205>.

</div>

<div id="ref-Gor25" class="csl-entry">

———. 2025. “Robust Estimation for Threshold Autoregressive
Moving-Average Models.” *Journal of Business and Economic Statistics*.
(.): in press. <https://doi.org/10.1080/07350015.2024.2412011>.

</div>

<div id="ref-Gor21" class="csl-entry">

Goracci, G., S. Giannerini, K.-S. Chan, and H. Tong. 2021. “Testing for
Threshold Effects in the TARMA Framework.” University of Bologna, Free
University of Bolzano, University of Iowa, London School of Economics.
<https://arxiv.org/abs/2103.13977>.

</div>

<div id="ref-Gor23" class="csl-entry">

———. 2023. “Testing for Threshold Effects in the TARMA Framework.”
*Statistica Sinica* 33 (3): 1879–1901.
https://doi.org/<https://doi.org/10.5705/ss.202021.0120>.

</div>

</div>
