# getPCA.test is used to compute the principal component scores
# for a given test sample
getPCA_test <- function(object, data){

  # object        : getPCA object
  # data          : functional data

  bs_basis <- object$bs_basis
  PCAcoef <- object$PCAcoef
  gp <- object$gp
  method.type <- object$method.type
  mean.tr <- c(object$meanScore$coefs)

  n <- dim(data)[1]
  p <- dim(data)[2]
  dimnames(data) = list(as.character(1:n), as.character(1:p))
  pcaobj <- smooth.basisPar(gp, t(data), bs_basis, Lfdobj=NULL, lambda=0)$fd

  if(method.type == "classical"){
    sdata <- scale(t(pcaobj$coefs), center = mean.tr, scale = FALSE)
    pcaobj2 <- pcaobj
    pcaobj2$coefs <- t(sdata)
    PCAscore_test = inprod(pcaobj2, PCAcoef)
    colnames(PCAscore_test) = 1:dim(PCAcoef$coefs)[2]
    for(i in 1:dim(PCAcoef$coefs)[2])
      colnames(PCAscore_test)[i] = paste("Score", i, sep = "")
  }else if(method.type == "robust"){
    sdata <- scale(t(pcaobj$coefs), center = mean.tr, scale = FALSE)
    pcaobj2 <- pcaobj
    pcaobj2$coefs <- t(sdata)
    PCAscore_test = inprod(pcaobj2, PCAcoef)
    colnames(PCAscore_test) = 1:dim(PCAcoef$coefs)[2]
    for(i in 1:dim(PCAcoef$coefs)[2])
      colnames(PCAscore_test)[i] = paste("Score", i, sep = "")
  }

  return(PCAscore.test = PCAscore_test)
}
