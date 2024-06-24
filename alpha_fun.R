# alpha_fun is used to compute the weight matrix in the M-estimator
alpha_fun <- function(z, c){
  # z       : standardized residuals
  # c       : tuning parameter
  
  ifelse(z == 0, 1, psi_fun(z, c) * z^{-1})
}