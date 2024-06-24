\name{fsac_pca}
\alias{fsac_pca}
\title{Scalar-on-function linear quantile regression based on functional principal component analysis}
\description{This function can be used to perform spatial scalar-on-function linear regression model \deqn{
Y = \rho Y W + \int X(t) \beta(t) dt + \epsilon} based on both the usual and robust functional principal component decomposition of the functional predictor.
}
\usage{
fsac_pca(y, x, nbasis = NULL, gpx = NULL, wei_mat,
         method.type = c("classical", "robust"),
         c1=1.4, c2=2.4, c3=1.65, tol = 1e-04)
}
\arguments{
  \item{y}{An \eqn{n \times 1}-dimensional matrix containing the observations of scalar response \eqn{Y}, where \eqn{n} denotes the sample size.}
  \item{x}{A matrix with dimension \eqn{n \times j} containing the observations of functional predictor \eqn{X(t)}, where \eqn{n} is the sample size and \eqn{j} denotes the number of grid points for \eqn{X(t)}.}
  \item{nbasis}{A numeric value denoting the number of B-spline basis expansion functions to be used to approximate the functional principal components.}
  \item{gpx}{A vector containing the grid points of the functional predictor \eqn{X(t)}.}
  \item{wei_mat}{An \eqn{n \times n} dimensional spatial weight matrix.}
  \item{method.type}{Method type used to estimate the spatial scalar-on-function linear regression model. Possibilities are "classical" and "robust".}
  \item{c1}{Tuning parameter.}
  \item{c2}{Tuning parameter.}
  \item{c3}{Tuning parameter.}
  \item{tol}{A numeric value used for the precision of the coefficient estimate.}
}
\value{A list object with the following components:
\item{b}{A vector containing the estimate of parameters of the regression model coducted between the scalar response and principal component scores of the functional predictor.}
\item{b0}{A vector containing the estimate of intercept parameter of the regression model coducted between the scalar response and principal component scores of the functional predictor.}
\item{bhat}{A vector containing the estimate of regression coefficient function.}
\item{rho}{Estimated spatial autocorrelation parameter.}
\item{sig}{Estimated standar deviation for the model.}
\item{fitted.values}{An \eqn{n \times 1}-dimensional matrix containing the fitted values of the scalar response.}
\item{residuals}{An \eqn{n \times 1}-dimensional matrix containing the residuals.}
\item{fpca}{A list object containing model details, such as number of basis functions, number of principal components, and grid points used for the functional predictor variable.}
\item{ncomp}{An integer denoting the number of principal components used in the analysis.}
}
\author{Ufuk Beyaztas, Abhijit Mandal and Han Lin Shang}
\examples{
# Generate data
sim.dat <- data_generation(n = 100, j = 101, rho = 0.5, sig.e = 1, out.p = 0.05)
y <- sim.dat$y
x <- sim.dat$x
wei_mat <- sim.dat$w

# Spatial functional regression model based on FPCA decomposition
fsacPCA <- fsac_pca(y=y, x=x, wei_mat = wei_mat, method.type = "classical")
# Robust spatial functional regression model based on FPCA decomposition
RfsacPCA <- fsac_pca(y=y, x=x, wei_mat = wei_mat, method.type = "robust")
}
