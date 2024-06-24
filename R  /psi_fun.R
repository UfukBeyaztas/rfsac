# Huber function
psi_fun <- function(z, c){
  # z       : standardized residuals
  # c       : tuning parameter
  
  z * pmin(1, c / abs(z))
}
