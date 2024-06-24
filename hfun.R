# h tilde function
hfun <- function(c){
  # c       : tuning parameter
  
  2 * c^2 * (1 - pnorm(c)) - 2 * c * dnorm(c) - 1 + 2 * pnorm(c)
}