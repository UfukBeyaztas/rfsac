\name{predict_flm_pca}
\alias{predict_flm_pca}
\title{Prediction for a scalar-on-function linear regression model based on functional principal component analysis}
\description{This function is used to make prediction for a new set of functional predictors based upon a fitted scalar-on-function linear regression model in the output of \code{\link{flm_pca}}.}
\usage{
predict_flm_pca(object, xnew)
}
\arguments{
  \item{object}{An output object obtained from \code{\link{flm_pca}}.}
  \item{xnew}{A matrix consisting of the new observations of functional predictor. The argument \code{xnew} must have the same length and the same structure as the input \code{x} of \code{\link{flm_pca}}.}
}
\value{An \eqn{n_{test} \times 1}-dimensional matrix of predicted values of the scalar response variable for the given set of new functional predictor \code{xnew}. Here, \eqn{n_{test}}, the number of rows of the matrix of predicted values, equals to the number of rows of \code{xnew}.
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

flmPCA <- flm_pca(y=y, x=x)
yhat_flm <- predict_flm_pca(object = flmPCA, xnew = x.test)
}
