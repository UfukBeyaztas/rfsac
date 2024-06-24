# cov_fun is used to compute the covariance between the response
# and predictor
cov_fun <- function(y, x){
  
  # y         : response
  # x         : predictor
  
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  cf <- numeric()
  for(i in 1:p){
    mod.i <- lm(y~x[,i])
    cf[i] <- mod.i$coefficients[2]
  }
  
  output <- diag(cov(x)) * cf
  
  return(output)
}