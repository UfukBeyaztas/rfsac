# fsac_pca is used to compute the maximum likelihood estimates or M-estimates 
# of the functional spatial AR model based on functional PCA decomposition.
fsac_pca <- function(y, x, nbasis = NULL, gpx = NULL, wei_mat,
                     method.type = c("classical", "robust"),
                     c1=1.4, c2=2.4, c3=1.65, tol = 1e-04){
  
  # y         : scalar response
  # x         : a matrix containing functional predictor
  # nbasis    : a vector containing number of basis functions
  # gpx       : a vector containing grid points 
  # wei_mat   : spatial weight matrix
  # method.type: type of estimation model; classical or robust
  # c1, c2, c3: numeric values used in the psi function
  # tol       : tolerance used in the iterative process
  
  method.type <- match.arg(method.type)
  n <- length(c(y))
  px <- dim(x)[2]
  
  if(is.null(gpx))
    gpx <- seq(1/px, 1-1/px, len = px)
  
  if(is.null(nbasis))
    nbasis <- round(min(n/4, 40))
  
  fpca <- getPCA(x, nbasis, gpx, method.type)
  fcomp <- fpca$PCAcoef
  fsco <- fpca$PCAscore
  evalbase <- fpca$evalbase
  
  if(method.type == "classical"){
    sac_model <- sac_ml(y, fsco, wei_mat)
    b0 <- sac_model$b0
    b <- sac_model$b
    rho <- sac_model$rho
    w <- sac_model$w
    sig <- sqrt(sac_model$sig2)
  }else if(method.type == "robust"){
    # Initial model and the parameters
    init.model <- sac_ml(y, fsco, wei_mat)
    x.k <- cbind(1, fsco)
    b0 <- init.model$b0
    b <- init.model$b
    rho <- init.model$rho
    sig <- sqrt(init.model$sig2)
    b. <-c(b0, b)
    
    # Iterative estimation algorithm
    
    repeat{
      
      # beta
      z <- as.vector((diag(n) - rho * wei_mat) %*% y - x.k %*% b.) / sig
      D <- diag(alpha_fun(z, c1))
      br <- ginv(t(x.k) %*% D %*% x.k) %*% t(x.k) %*% D %*% (y - rho * wei_mat %*% y)
      
      # sigma
      zr1 <- ((diag(n) - rho * wei_mat) %*% y - x.k %*% br) / sig
      sigr <- sig^2 / (n * hfun(c2)) * 
        t(as.matrix(psi_fun(zr1, c2))) %*% as.matrix(psi_fun(zr1, c2))
      sigr <- as.numeric(sqrt(sigr))
      
      # rho
      temp <- function(u){
        zrho <- ((diag(n) - u * wei_mat) %*% y - x.k %*% br) / sigr
        Grho <- wei_mat %*% ginv(diag(n) - u * wei_mat)
        buff <- (1/sigr) * t(Grho %*% x.k %*% br) %*% psi_fun(zrho, c3) +
          t(psi_fun(zrho, c3)) %*% t(Grho) %*% psi_fun(zrho, c3) - 
          sum(diag(Grho)) * hfun(c3)
        return(buff^2)
      }
      
      
      # set the lower bound for rho
      # eigsW <- eigen(wei_mat)$values
      # lb <- -abs(1/min(eigsW))
      # if(lb < -.99) lb <- -.99
      # For general use, we set lb <- -.99
      
      rhor <- optimize(temp, lower=-.99, upper=.99)$minimum
      
      # if the lower bound determined by the minimum of eigenvalue use following
      # rhor <- optimize(temp, lower=lb, upper=.99)$minimum
      
      b.diff <- max(abs(br-b.))
      sig.diff = abs(sigr - sig)
      rho.diff = abs(rhor - rho)
      
      if(b.diff < tol || sig.diff < tol || rho.diff < tol)
        break
      
      b. <- br
      sig <- sigr
      rho <- rhor
    }
    
    b0 <- b.[1]
    b <- b.[-1]
  }
  
  fits <- ginv(diag(n) - rho*wei_mat) %*% as.matrix(rep(1,n)) * b0 + 
    ginv(diag(n) - rho*wei_mat) %*% fsco %*% b
  
  b_hat_t <- evalbase %*% (fcomp$coefs %*% b)
  resids <- y - fits
  
  return(list(b = b, b0 = b0,
              bhat = b_hat_t,
              rho = rho,
              sig = sig,
              fitted.values = fits,
              residuals = resids,
              fpca = fpca,
              ncomp = fpca$ncomp))
}
