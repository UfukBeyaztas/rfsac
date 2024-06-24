\name{fsac_pls}
\alias{fsac_pls}
\title{Scalar-on-function linear quantile regression based on functional partial least squares}
\description{This function can be used to perform spatial scalar-on-function linear regression model \deqn{
Y = \rho Y W + \int X(t) \beta(t) dt + \epsilon} based on both the usual and robust functional partial least squares decomposition of the functional predictor.
}
\usage{
fsac_pls(y, x, h, nbasis = NULL, gpx = NULL, wei_mat,
         method.type = c("classical","robust"),
         probp1 = 0.95, hampelp2 = 0.975, hampelp3 = 0.999,
         maxit = 1000, conv = 0.01,
         c1=1.4, c2=2.4, c3=1.65, tol = 1e-04)
}
\arguments{
  \item{y}{An \eqn{n \times 1}-dimensional matrix containing the observations of scalar response \eqn{Y}, where \eqn{n} denotes the sample size.}
  \item{x}{A matrix with dimension \eqn{n \times j} containing the observations of functional predictor \eqn{X(t)}, where \eqn{n} is the sample size and \eqn{j} denotes the number of grid points for \eqn{X(t)}.}
  \item{h}{A numeric value denoting the number of functional partial least regression components to be computed.}
  \item{nbasis}{A numeric value denoting the number of B-spline basis expansion functions to be used to approximate the functional principal components.}
  \item{gpx}{A vector containing the grid points of the functional predictor \eqn{X(t)}.}
  \item{wei_mat}{An \eqn{n \times n} dimensional spatial weight matrix.}
  \item{method.type}{Method type used to estimate the spatial scalar-on-function linear regression model. Possibilities are "classical" and "robust".}
  \item{probp1}{A numeric value used to determine the first outlier cutoff point for the weights.}
  \item{hampelp2}{A numeric value used to determine the first outlier cutoff point for the weights.}
  \item{hampelp3}{A numeric value used to determine the third outlier cutoff point for the weights.}
  \item{maxit}{An integer value defining the maximum iteration used to achieve convergence.}
  \item{conv}{A numeric value used for the precision of the coefficient estimate.}
  \item{c1}{Tuning parameter.}
  \item{c2}{Tuning parameter.}
  \item{c3}{Tuning parameter.}
  \item{tol}{Tolerance for the approximation of M-estimation.}
}
\value{A list object with the following components:
\item{b}{A vector containing the estimate of parameters of the regression model coducted between the scalar response and partial least squares scores of the functional predictor.}
\item{b0}{A vector containing the estimate of intercept parameter of the regression model coducted between the scalar response and partial least squares scores of the functional predictor.}
\item{bhat}{A vector containing the estimate of regression coefficient function.}
\item{rho}{Estimated spatial autocorrelation parameter.}
\item{sig}{Estimated standar deviation for the model.}
\item{fitted.values}{An \eqn{n \times 1}-dimensional matrix containing the fitted values of the scalar response.}
\item{residuals}{An \eqn{n \times 1}-dimensional matrix containing the residuals.}
\item{fdd}{A list object containing model details, such as number of basis functions, number of partial least squares components, and grid points used for the functional predictor variable.}
}
\author{Ufuk Beyaztas, Abhijit Mandal and Han Lin Shang}
\examples{
# Generate data
sim.dat <- data_generation(n = 100, j = 101, rho = 0.5, sig.e = 1, out.p = 0.05)
y <- sim.dat$y
x <- sim.dat$x
wei_mat <- sim.dat$w

fsacPCA <- fsac_pca(y=y, x=x, wei_mat = wei_mat, method.type = "classical")

# Spatial functional regression model based on FPLS decomposition
fsacPLS <- fsac_pls(y=y, x=x, h = fsacPCA$ncomp, wei_mat = wei_mat, method.type = "classical")
# Robust spatial functional regression model based on FPLS decomposition
RfsacPLS <- fsac_pls(y=y, x=x, h = fsacPCA$ncomp, wei_mat = wei_mat, method.type = "robust")
}
