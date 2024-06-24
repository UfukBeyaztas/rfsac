\name{predict_fsac_pls}
\alias{predict_fsac_pls}
\title{Prediction for a spatial scalar-on-function linear regression model based on functional partial least squares}
\description{This function is used to make prediction for a new set of spatial weight matrix and functional predictors based upon a fitted spatial scalar-on-function linear regression model in the output of \code{\link{fsac_pls}}.
}
\usage{
predict_fsac_pls(object, xnew, wnew)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{fsac_pls}}.}
  \item{xnew}{A matrix consisting of the new observations of functional predictor. The argument \code{xnew} must have the same length and the same structure as the input \code{x} of \code{\link{fsac_pls}}.}
  \item{wnew}{A matrix consisting of spatial weight matrix for the test set.}
}
\value{An \eqn{n_{test} \times 1}-dimensional matrix of predicted values of the scalar response variable for the given spatial weight matrix and a set of new functional predictor \code{xnew}. Here, \eqn{n_{test}}, the number of rows of the matrix of predicted values, equals to the number of rows of \code{xnew}.
}
\author{Ufuk Beyaztas, Abhijit Mandal and Han Lin Shang}
\examples{
# Generate data
sim.dat <- data_generation(n = 100, j = 101, rho = 0.5, sig.e = 1, out.p = 0.05)
y <- sim.dat$y
x <- sim.dat$x
wei_mat <- sim.dat$w

# Generate test data
test.dat <- data_generation(n = 100, j = 101, rho = 0.5, sig.e = 1)
x.test <- test.dat$x
y.test <- test.dat$y
w.test <- test.dat$w

fsacPCA <- fsac_pca(y=y, x=x, wei_mat = wei_mat, method.type = "classical")

# Spatial functional regression model based on FPCA decomposition
fsacPLS <- fsac_pls(y=y, x=x, h = fsacPCA$ncomp, wei_mat = wei_mat, method.type = "classical")
yhat_PLS <- predict_fsac_pls(object = fsacPLS, xnew = x.test, wnew = w.test)
# Robust spatial functional regression model based on FPCA decomposition
RfsacPLS <- fsac_pls(y=y, x=x, h = fsacPCA$ncomp, wei_mat = wei_mat, method.type = "robust")
yhat_RPLS <- predict_fsac_pls(object = RfsacPLS, xnew = x.test, wnew = w.test)
}
