## Declare ACValues as global to avoid CRAN check note

globalVariables(c("ACValues"))

## Given pa and p, gives the corresponding tabulated critical values
## for the ARMA vs TARMA test
## based on ACValues
CritValues <- function(pa,p){
    ind <- which.min(abs(ACValues[,'pi']-pa))  # selects the row such that pi is the closest to pa
    ACValues[ind,paste(p,c(90,95,99),sep='-')] # selects the critical values based on p
}

