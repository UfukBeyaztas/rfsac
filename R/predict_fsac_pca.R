# predict_fsac_pca is used to obtain the predictions for a given
# new functional predictor and spatial weight matrix based on fsac_pca model
predict_fsac_pca <- function(object, xnew, wnew){

  # object    : fsac_pca object
  # xnew      : a matrix containing functional predictor in the test sample
  # wnew      : spatial weight matrix for the test sample

  n <- dim(xnew)[1]
  rho <- object$rho
  b0 <- object$b0
  b <- object$b
  fpca <- object$fpca
  method.type <- fpca$method.type
  fsco.test <- getPCA_test(fpca, xnew)

  predicted.values <- ginv(diag(n) - rho*wnew) %*% as.matrix(rep(1,n)) * b0 +
    ginv(diag(n) - rho*wnew) %*% fsco.test %*% b

  return(predicted.values)
}
