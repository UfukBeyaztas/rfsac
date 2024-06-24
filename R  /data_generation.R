data_generation = function(n, j, rho, sig.e, out.p = 0){
  
  s = seq(0, 1, length.out = j)
  
  ksi = list()
  for(ik in 1:5){
    ksi[[ik]] = rnorm(n, 0, sd = (4*ik^(-3/2)))
  }
  
  phi = list()
  for(ik in 1:5){
    phi[[ik]] = sin(ik * pi * s) - cos(ik * pi * s)
  }
  
  fX = Reduce("+", lapply(1:5, function(k){ksi[[k]] %*% t(phi[[k]])}))
  
  vBeta = sin(2*pi * s)
  
  fX = fdata(fX, argvals = s)
  vBeta = fdata(vBeta, argvals = s)
  
  err = rnorm(n, mean=0, sd=sig.e)
  
  argx = inprod.fdata(fX, vBeta)
  
  wei <- matrix(, n, n)
  
  for(i in 1:n){
    for(j in 1:n){
      if(i != j){
        wei[i,j] <- 1/(abs(i-j))
      }else{
        wei[i,j] <- 0
      }
    }
  }
  
  W <- matrix(0, n, n)
  for(i in 1:n)
    W[i,] <- wei[i,] / sum(wei[i,])
  
  
  fYe = ginv(diag(n) - rho*W) %*% argx + ginv(diag(n) - rho*W) %*% err
  out.index <- NULL
  
  if(out.p > 0){
    
    nout <- round(n * out.p)
    out.index <- sample(1:n, nout)
    err[out.index] <- rnorm(nout, mean=10, sd=1)
    
    ksi.out = list()
    for(ik in 1:5){
      ksi.out[[ik]] = rnorm(n, 1, sd = (4*ik^(-3/2)))
    }
    
    phi.out = list()
    for(ik in 1:5){
      phi.out[[ik]] = sin(ik * pi * s) - cos(ik * pi * s)
    }
    
    fX.out = Reduce("+", lapply(1:5, function(k){ksi.out[[k]] %*% t(phi.out[[k]])}))
    
    fX1 <- fX$data
    fX1[out.index,] <- fX.out[out.index,]
    fX = fdata(fX1, argvals = s)
    argx = inprod.fdata(fX, vBeta)
    fYe = ginv(diag(n) - rho*W) %*% argx + ginv(diag(n) - rho*W) %*% err
  }
  
  return(list("y" = fYe, "x" = fX$data, w = W, tcoefs = vBeta, out.index = out.index))
  
}
