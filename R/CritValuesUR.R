## Declare supLMQur as global to avoid CRAN check note

globalVariables("supLMQur")

CritValuesUR <- function(n,pa,th){
  ## INPUT 
  ## n : sample size
  ## pa: lower threshold level, 
  ## th: MA(1) parameter , 
  ## OUTPUT 
  ## tabulated critical values for the supLM Unit Root test IMA vs TARMA 
  ## based on supLMQur
  pset    <- c(0.01,seq(from=0.05,to=0.40,by=0.05))
  npar    <- length(pset)
  lserie  <- c(seq(100,1000,by=100),5000)
  nset    <- length(lserie)
  thset   <- c(-0.9,-0.6,-0.3,0,0.3,0.6,0.9)
  
  i.pa <- which.min(abs(pset-pa))  # the closest to pa
  i.th <- which.min(abs(thset-th)) # the closest to th
  i.n  <- which.min(abs(lserie-n)) # the closest to n
  return(supLMQur[,i.pa,i.th,i.n])
}