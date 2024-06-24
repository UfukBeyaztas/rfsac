# predict_flm_pca is used to obtain the predictions for a given
# new functional predictor based on flm_pca model
predict_flm_pca <- function(object, xnew){

  # object    : flm_pca object
  # xnew      : a matrix containing functional predictor in the test sample

  n <- dim(xnew)[1]
  b0 <- object$b0
  b <- object$b
  fpca <- object$fpca
  method.type <- fpca$method.type
  fsco.test <- getPCA_test(fpca, xnew)

  predicted.values <- b0 + fsco.test %*% b

  return(predicted.values)
}
