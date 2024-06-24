# sac_ml is used to compute the maximum likelihood estimates
# of the spatial AR model.
sac_ml <- function(y, x, wei_mat){
  
  # y             : response variable
  # x             : predictor matrix
  # wei_mat       : spatial weight matrix
  
  W <- as(wei_mat,"dgCMatrix")
  
  # compute the vector of eigenvalues
  eigvar <- eigen(W, only.values=TRUE, symmetric=TRUE)$values
  
  wy <- as.vector(wei_mat %*% y)
  xmat <- cbind(1, x, wy)
  n <- dim(xmat)[1]
  p <- dim(xmat)[2]
  
  # concentrated likelihood
  xm_rho <- cbind(1, x)
  beta_rho <- ginv(t(xm_rho) %*% xm_rho) %*% t(xm_rho) %*% y
  err_rho <- y - cbind(1, x) %*% beta_rho
  
  beta_wy <- ginv(t(xm_rho) %*% xm_rho) %*% t(xm_rho) %*% wy
  err_wy <- wy - cbind(1, x) %*% beta_wy
  
  lc <- function(u){
    err_u <- (err_rho - u*err_wy)^2
    temp = n*log(mean(err_u))/2 - sum(log(1-u*eigvar))
    return(temp)
  }
  
  rho_init <- optimize(lc, lower=-.99, upper=.99)$minimum
  
  err <- y-rho_init*wy
  beta_init <- ginv(t(xm_rho) %*% xm_rho) %*% t(xm_rho) %*% err
  beta_hat <- c(beta_init, rho_init)
  sig2_init <- sum((err - cbind(1, x) %*% beta_init)^2) / (n - length(beta_init))
  
  # nlm to minimize log-likelihood function
  log_lik <- function(beta){
    s2 = beta[p+1]
    rho = beta[p]
    bmat <- beta[1:p]
    xb <- xmat%*%bmat
    e <- y - xb
    out <-  n*log(pi)/2 + n*log(s2)/2 + sum(e^2)/(2*s2) - sum(log(1 - rho*eigvar))
    return(out)
  }
  
  options(warn=-1)
  nlfit <- nlm(log_lik, c(beta_init, rho_init, sig2_init), hessian=TRUE)
  options(warn=0)
  
  best <- nlfit$estimate
  b0_hat <- best[1]
  b_hat <- best[2:(p-1)]
  rho_hat <- best[p]
  sig2_hat <- best[p+1]
  
  xtx <- crossprod(xmat)
  xtx[p, p] = xtx[p, p] + sig2_hat*sum((eigvar/(1-rho_hat*eigvar))^2)
  vmat <- sig2_hat*ginv(xtx)
  
  return(list(b0 = b0_hat,
              b = b_hat,
              rho = rho_hat,
              sig2 = sig2_hat,
              w = wei_mat,
              vmar = vmat))
}