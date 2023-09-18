
# tseriesTARMA changelog

## 0.3-4

- Removed a leftover call to the `Fortran` intrinsic `RANDOM_NUMBER` function.

- Fixed the docs to cope with the roxygen2 issue https://github.com/r-lib/roxygen2/issues/1491


## 0.3-3

- Now the residuals in `TARMA.fit` inherit the time attributes of the original series.

- Fixed `print.TARMAtest` when `fixed` is used to fix parameter values in `TARMA.fit2`.

- In `TARMAGARCH.test`: set `include.mean=TRUE` in the ARMA fit and removed centering.

- Minor fixes to the documentation of `TARMAGARCH.test`.

- Updated the e-mail address of Greta Goracci in the docs.

- Updated the Statistica Sinica 2023 reference.

## 0.3-2

- Fixed missing dependency in MAKEVARS: the Fortran library TARMAurtest depends on the Fortran module TARMA_MOD

- Fixed broken references in the `README.md` file. Now `README.Rmd` is used to generate it.

- Fixed missing `PACKAGE` in several `.Fortran` calls.


## 0.3-1

- First CRAN submission.

- Fixed urls and DESCRIPTION file to avoid notes from CRAN. 

- Added `value` field to `plot.tsfit`, `predict.TARMA`, `print.TARMA`, `print.TARMAtest`. 

- Improved the examples of `predict.TARMA`.

- Removed `\dontrun` or replaced it with `\donttest`.
 
## 0.3-0

- The new functions `TARMAur.test` and `TARMAur.B.test` implement the asymptotic and bootstrap supLM unit root test for an integrated MA against a stationary TARMA process. The print method has been updated to handle them and uses tabulated critical values made available through the dataset `supLMQur`.  

- Added the Journal of Econometrics reference for `TAR.test` and `TAR.test.B`.

- The new function `TAR.test` implements the robust supLM test for TAR nonlinearity.

- The dataset `ACValues` contains the tabulated critical values of Andrews(2003), Econometrica. 

- `TAR.test`, `TAR.test.B`, `TARMA.test`, `TARMAGARCH.test` have been modified to use the new class `TARMAtest`. The print method for `TARMAtest` is implemented and extracts and prints the critical values from the dataset `ACValues`.

- `TARMA.test.B` has been replaced by the new function `TAR.test.B`. Added the bibtex item for the relevant paper.

- The Fortran module `TARMA_MOD` has been updated. Also, the `ISO_C_BINDING` facilities for floating point types are used.  
  
- Removed the switch `shift` from `TARMA.test` and `TARMAGARCH.test`.

- `TARMAGARCH.test` now always uses `method` for arima fitting. 

- Minor cosmetic polishing of the documentation of `TARMA.test` and `TARMAGARCH.test`.

- The routines `tarma.fit` and `tarma.fit2` have been renamed`TARMA.fit` and `TARMA.fit2`, for consistency.

## 0.2-5

- In `tarma2.fit` a warning is issued if the regressors and the series have different lengths.

- Fixed a bug in `tarma2.fit`. The residuals lost the `tsp` attributes of the original series.

- Fixed a bug in `tarma.fit`. Now the robust standard errors routine for the `robust` method works also for TMA order greater than one. The internal routine D2eps.R has been fixed.

- The documentation of `tarma.fit` has been improved. Most notably, the examples section has been expanded and is now run. 

- Fixed the examples of `plot.tsfit` and `predict.TARMA`: they had not been updated after the change to `tarma.fit` in version `0.1-4`.

- `print.TARMA` now reports also the value of `alpha` and `qu` if the estimation method is `robust`.

- `plot.tsfit` gains more flexibility by allowing control over: line type with `ltype`, line width with `lwdt`.
Also, the default line colours and widths have been changed. The examples are now run.

- Minor improvements to the docs: the defaults of `predict.TARMA` have been changed. Added a reference for the robustness paper.   

## 0.2-4

- Now `tarma.fit2` allows to specify different subsets of lags for the two regimes: `tar1.lags` and `tar2.lags`. `tarma.fit` has been changed accordingly and now passes the exact subsets to `tarma.fit2` for initial estimation. The documentation of `tarma.fit2` has been updated.

- `print.TARMA` has been improved and the methods for the TARMA class have been moved to a separate file.

- `plot.tsfit` now uses `inherits` to avoid the note from checking.

- Minor corrections and improvements to the documentation of `tarma.fit`.

- Minor improvements of `Description` and `Readme`.

## 0.2-3

- Fixed a bug in `predict.TARMA`: the number of discarded observations did not take into account `ar.lags`. 

- `TARMA.test`  and `TARMA.test.B` now always use `method` for arima fitting. 

- `TARMA.test`  and `TARMA.test.B` now also work with `zoo` objects. 

## 0.2-2

- The robust version of `tarma.fit` has been polished and is no more experimental. 
    Fixed a bug in the derivation of robust standard errors.

- Fixed a in bug in `tarma.fit` for which the trimmed estimation that was meant to be used only for the initial estimation was used also in the subsequent steps. Fixed alpha in line 407.

## 0.2-1

- The robust version of `tarma.fit` has been completely rewritten and uses an experimental additional iteratively reweighted least squares step to estimate the robust weights. 

## 0.2-0 

- The routines `ARMAvsTARMA.test`  and `ARMAvsTARMA.B.test` have ben renamed  `TARMA.test`  and `TARMA.B.test` , respectively.

- Added the routine `TARMAGARCH.test` .

## 0.1-5:

- The method `print.TARMA` has been rewritten.

- The docs have been rewritten to use the package `mathjaxr`.

- Fixed a bug in line 259 of `tarma.fit2`: replaced `>` with `!=`.

- The output slot `method` is added to `tarma.fit2` for consistency with `tarma.fit`. It is set 
    to `MLE` and has no effect other than allowing to differentiate between the two routines.

- The output slots `arlag` and `include.int` are added to `tarma.fit2`.

## 0.1-4:

- `tarma.fit` has been rewritten. Now it can fit a full subset TARMA(p1,p2,q1,q2) model.
    Some ergodicity bounds for the optimization routines have been also added. 
    The objective functions have been rewritten in a recursive fashion and the code now is 
    several times faster than before. `predict.TARMA` has been modified to cope with the new routine.

- Added a warning upon failed convergence in `tarma.fit`.

- Added the S3 method `residuals.TARMA`.

- `tarma.fit` now outputs the values of the AIC and BIC computed 
    according to Chan & Cryer (2008) and Chen & Tsay (2019).

- The output of `tarma.fit2` has been modified to match that of the new version of `tarma.fit`.

- `ARMAvsTARMA.test` now includes the full test for ARMA versus TARMA when `ma.fixed=FALSE`. 
   The output of the test now includes the information on the parameters tested. The documentation
   has been amended accordingly.
   
- Fixed a bug in `predict.TARMA`: now it handles correctly time series objects 
    with `frequency` different from 1.

- Now `tarma.sim` simulates from a full TARMA(p1,p2,q1,q2) model where q1 and q2 can be different.


## 0.1-3:

- `tarma.fit2` now gains the possibility to add a seasonal MA part. The output has been amended accordingly.

- Fixed a bug in `tarma.sim` when the (T)MA order was greater than the (T)AR order.

- Now the output of `tarma.fit2` has no `NA`s and can be used in `tarma.sim`.

## 0.1-2:

- Changed the option `fitted` in `plot.fit` in `plot.tsfit`.

- Fixed a bug in `tarma.sim`: now it handles correctly time series objects with `frequency` different from 1.

- `ARMAvsTARMA.test` now includes the test for AR versus TAR when `ma.ord=0`. 

- Fixed a bug in the output of `ARMAvsTARMA.test.B`: `Tb` was outside the output list.

## 0.1-1:

- `tarma.fit2` now allows the inclusion of covariates in the model.

## 0.1-0:

- First version
