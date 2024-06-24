# predict_fsac_pls is used to obtain prediction for a given new functional predictor
# and spatial weight matrix based on fsac_pls model

predict_fsac_pls <- function(object, xnew, wnew){

  # object    : fsac_pls object
  # xnew      : a matrix containing functional predictor in the test sample
  # wnew      : spatial weight matrix for the test sample

  n <- dim(xnew)[1]
  rho <- object$rho
  b0 <- object$b0
  nbasis <- object$fdd$nbasis
  gpx <- object$fdd$gpx
  m.tr <- object$fdd$m.tr
  fin.cf <- object$fdd$fin.cf

  BS.sol.test <- getAmat(data = xnew, nbasis = nbasis, gp = gpx)
  Amat.test <- BS.sol.test$Amat
  Amat.test <- scale(Amat.test, center = m.tr, scale = F)

  predicted.values <- ginv(diag(n) - rho*wnew) %*% as.matrix(rep(1,n)) * b0 +
    ginv(diag(n) - rho*wnew) %*% Amat.test %*% fin.cf

  return(predicted.values)
}
