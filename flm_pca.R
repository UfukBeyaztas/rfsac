# flm_pca is used to compute the maximum likelihood estimates for  
# functional linear regression based on the functional PCA decomposition.
flm_pca <- function(y, x, nbasis = NULL, gpx = NULL){
  
  # y         : scalar response
  # x         : a matrix containing functional predictor
  # nbasis    : a vector containing number of basis functions
  # gpx       : a vector containing grid points 
  
  n <- length(c(y))
  px <- dim(x)[2]
  
  if(is.null(gpx))
    gpx <- seq(1/px, 1-1/px, len = px)
  
  if(is.null(nbasis))
    nbasis <- round(min(n/4, 40))
  
  fpca <- getPCA(x, nbasis, gpx, method.type = "classical")
  fcomp <- fpca$PCAcoef
  fsco <- fpca$PCAscore
  evalbase <- fpca$evalbase
  
  
  xmk <- cbind(1, fsco)
  cfs <- ginv(t(xmk) %*% xmk) %*% t(xmk) %*% y
  b0 <- cfs[1]
  b <- cfs[-1]
  fits <- cbind(1, fsco) %*% cfs
  
  b_hat_t <- evalbase %*% (fcomp$coefs %*% b)
  resids <- y - fits
  sig <- sqrt(as.numeric((t(resids) %*% resids)/(n-2)))
  
  return(list(b = b, b0 = b0,
              bhat = b_hat_t,
              fitted.values = fits,
              residuals = resids,
              sig = sig,
              fpca = fpca))
}