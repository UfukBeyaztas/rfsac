# fsac_pls is used to compute the maximum likelihood estimates
# of the functional spatial AR model based on functional PLS decomposition
fsac_pls <- function(y, x, h, nbasis = NULL, gpx = NULL, wei_mat,
                     method.type = c("classical","robust"),
                     probp1 = 0.95, hampelp2 = 0.975, hampelp3 = 0.999,
                     maxit = 1000, conv = 0.01,
                     c1=1.4, c2=2.4, c3=1.65, tol = 1e-04){
  
  # y         : scalar response
  # x         : a vector containing functional predictor
  # h         : number of iterations
  # nbasis    : number of basis functions
  # gpx       : a vector containing grid points 
  # wei_mat   : spatial weight matrix
  
  method.type <- match.arg(method.type)
  
  n <- length(c(y))
  px <- dim(x)[2]
  
  if(is.null(gpx))
    gpx <- seq(1/px, 1-1/px, len = px)
  
  if(is.null(nbasis))
    nbasis <- round(min(n/4, 40))
  
  
  x0 <- x
  y0 <- y
  
  BS.sol <- getAmat(data = x0, nbasis = nbasis, gp = gpx)
  x <- BS.sol$Amat
  evalbase <- BS.sol$evalbase
  sinp_mat <- BS.sol$sinp_mat
  
  yx.s <- scale_fun(data = cbind(y, x), method.type = method.type)
  xs <- yx.s[,-1]
  ys <- as.matrix(yx.s[,1])
  
  if(method.type == "classical"){
    m.pls <- pls_fun(y = ys, x = xs, h = h)
  }else if(method.type == "robust"){
    
    wx <- sqrt(apply(yx.s[,-1]^2, 1, sum))
    wx <- wx / median(wx)
    wy <- abs(yx.s[,1])
    zerows <- vector(length=0)
    
    if(length(wy) / 2 > sum(wy == 0)){
      wy <- wy / (1.4826 * median(wy))
    }else{
      wy <- wy / (1.4826 * median(wy[wy != 0]))
    }
    
    probct <- qnorm(probp1)
    hampelb <- qnorm(hampelp2)
    hampelr <- qnorm(hampelp3)
    wx[which(wx <= probct)] <- 1
    wx[which(wx > probct & wx <= hampelb)] <- probct / abs(wx[which(wx > probct & wx <= hampelb)])
    wx[which(wx > hampelb & wx <= hampelr)] <- probct * (hampelr - abs(wx[which(wx > hampelb & wx <= hampelr)]))/
      (hampelr - hampelb) * 1 / abs(wx[which(wx > hampelb & wx <= hampelr)])
    wx[which(wx > hampelr)] <- 0
    wy[which(wy <= probct)] <- 1
    wy[which(wy > probct & wy <= hampelb)] <- probct/abs(wy[which(wy > probct & wy <= hampelb)])
    wy[which(wy > hampelb & wy <= hampelr)] <- probct * (hampelr - abs(wy[which(wy > hampelb & wy <= hampelr)]))/
      (hampelr - hampelb) * 1 / abs(wy[which(wy > hampelb & wy <= hampelr)])
    wy[which(wy > hampelr)] <- 0
    
    w <- wx * wy
    
    if(any(w < 1e-6)){
      w0 <- which(w < 1e-6)
      w <- replace(w, list=w0, values = 1e-6)
      we <- w
    } else {
      wxe <- wx
      wye <- wy
      we <- w
    }
    
    yx.w <- as.data.frame(yx.s) * sqrt(we)
    loops <- 1
    rold <- 10^-5
    difference <- 1
    
    
    while((difference > conv) && loops < maxit){
      m.pls <- pls_fun(y = as.matrix(yx.w[,1]), x = as.matrix(yx.w[,-1]), h = h)
      yp <- m.pls$fitted.values
      r <- yx.s[,1] - yp
      b <- m.pls$d.coef
      Tpls <- as.data.frame(m.pls$T) / sqrt(we)
      
      if (length(r) / 2 > sum(r == 0)){
        r <- abs(r) / (1.4826 * median(abs(r)))
      }else{
        r <- abs(r) / (1.4826 * median(abs(r[r != 0])))
      }
      
      dt <- scale_fun(data = Tpls, method.type = method.type)
      wtn <- sqrt(apply(dt^2, 1, sum))
      wtn <- wtn/median(wtn)
      
      probct <- qnorm(probp1)
      hampelb <- qnorm(hampelp2)
      hampelr <- qnorm(hampelp3)
      wye <- r
      wye[which(r <= probct)] <- 1
      wye[which(r > probct & r <= hampelb)] <- probct / abs(r[which(r > probct & r <= hampelb)])
      wye[which(r > hampelb & r <= hampelr)] <- probct * (hampelr-abs(r[which(r > hampelb & r <= hampelr)]))/
        (hampelr - hampelb) * 1 / abs(r[which(r > hampelb & r <= hampelr)])
      wye[which(r > hampelr)] <- 0
      wye <- as.numeric(wye)
      
      probct <- qchisq(probp1, h)
      hampelb <- qchisq(hampelp2, h)
      hampelr <- qchisq(hampelp3, h)
      wte <- wtn
      wte[which(wtn <= probct)] <- 1
      wte[which(wtn > probct & wtn <= hampelb)] <- probct / abs(wtn[which(wtn > probct & wtn <= hampelb)])
      wte[which(wtn > hampelb & wtn <= hampelr)] <- probct * (hampelr-abs(wtn[which(wtn > hampelb & wtn <= hampelr)])) /
        (hampelr - hampelb) * 1 / abs(wtn[which(wtn > hampelb & wtn <= hampelr)])
      wte[which(wtn > hampelr)] <- 0
      
      difference <- abs(sum(b^2) - rold)/rold
      rold <- sum(b^2)
      we <- wye * wte
      
      if(any(we < 1e-6)){
        w0 <- which(we < 1e-6)
        we <- replace(we, list = w0, values = 1e-6)
        zerows <- unique(c(zerows, as.numeric(names(w0))))
      }
      
      if(length(zerows) >= (n/2)){
        break
      }
      
      yx.w <- as.data.frame(yx.s) * sqrt(we)
      
      loops <- loops + 1
    }
    
    if (difference > conv){
      warning(paste("Convergence was not achieved. The scaled difference between norms of the coefficient vectors is ",
                    round(difference, digits=4)))
    }
    
    w <- we
    w[zerows] <- 0
    wt <- wte
    wt[zerows] <- 0
    wy <- wye
    wy[zerows] <- 0
  }
  
  P <- m.pls$P
  R <- m.pls$R
  W <- m.pls$W
  if(method.type == "classical"){
    T <- m.pls$T
  }else if(method.type == "robust"){
    T <- xs %*% R
  }
  
  
  if(method.type == "classical"){
    sac_model <- sac_ml(y0, T, wei_mat)
    b0 <- sac_model$b0
    b <- sac_model$b
    rho <- sac_model$rho
    w <- sac_model$w
    sig <- sqrt(sac_model$sig2)
  }else if(method.type == "robust"){
    # Initial model and the parameters
    init.model <- sac_ml(y0, T, wei_mat)
    x.k <- cbind(1, T)
    b0 <- init.model$b0
    b <- init.model$b
    rho <- init.model$rho
    sig <- sqrt(init.model$sig2)
    b. <-c(b0, b)
    
    # Iterative estimation algorithm
    
    repeat{
      
      # beta
      z <- as.vector((diag(n) - rho * wei_mat) %*% y0 - x.k %*% b.) / sig
      D <- diag(alpha_fun(z, c1))
      br <- ginv(t(x.k) %*% D %*% x.k) %*% t(x.k) %*% D %*% (y0 - rho * wei_mat %*% y0)
      
      # sigma
      zr1 <- ((diag(n) - rho * wei_mat) %*% y0 - x.k %*% br) / sig
      sigr <- sig^2 / (n * hfun(c2)) * 
        t(as.matrix(psi_fun(zr1, c2))) %*% as.matrix(psi_fun(zr1, c2))
      sigr <- as.numeric(sqrt(sigr))
      
      # rho
      temp <- function(u){
        zrho <- ((diag(n) - u * wei_mat) %*% y0 - x.k %*% br) / sigr
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
  
  fits <- ginv(diag(n) - rho*w) %*% as.matrix(rep(1,n)) * b0 + 
    ginv(diag(n) - rho*w) %*% T %*% b
  
  b_hat_t <- evalbase %*% (t(solve(sinp_mat)) %*% W %*% b)
  resids <- y0 - fits
  
  fd_details <- list()
  fd_details$nbasis <- nbasis
  fd_details$gpx <- gpx
  fd_details$fin.cf <- R %*% as.matrix(b)
  if(method.type == "classical"){
    fd_details$m.tr <- apply(x, 2, mean)
  }else if(method.type == "robust"){
    fd_details$m.tr <- l1median(x)
  }
  
  return(list(b = b, b0 = b0,
              bhat = b_hat_t,
              rho = rho,
              sig = sig,
              fitted.values = fits,
              residuals = resids,
              fdd = fd_details))
  
}