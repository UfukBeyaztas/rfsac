\name{flm_pca}
\alias{flm_pca}
\title{Scalar-on-function linear regression based on functional principal component analysis}
\description{This function can be used to perform scalar-on-function linear regression model \deqn{
Y = \int X(t) \beta(t) dt + \epsilon,} based on the functional principal component decomposition of the functional predictor.
}
\usage{
flm_pca(y, x, nbasis = NULL, gpx = NULL)
}
\arguments{
  \item{y}{An \eqn{n \times 1}-dimensional matrix containing the observations of scalar response \eqn{Y}, where \eqn{n} denotes the sample size.}
  \item{x}{A matrix with dimension \eqn{n \times j} containing the observations of functional predictor \eqn{X(t)}, where \eqn{n} is the sample size and \eqn{j} denotes the number of grid points for \eqn{X(t)}.}
  \item{nbasis}{A numeric value denoting the number of B-spline basis expansion functions to be used to approximate the functional principal components for the functional predictor \eqn{X(t)}.}
  \item{gpx}{A vector containing the grid points of the functional predictor \eqn{X(t)}.}
}
\value{A list object with the following components:
\item{b}{A vector containing the estimate of parameters of the regression model coducted between the scalar response and principal component scores of the functional predictor.}
\item{b0}{A vector containing the estimate of intercept parameter of the regression model coducted between the scalar response and principal component scores of the functional predictor.}
\item{bhat}{A vector containing the estimate of regression coefficient function.}
\item{fitted.values}{An \eqn{n \times 1}-dimensional matrix containing the fitted values of the scalar response.}
\item{residuals}{An \eqn{n \times 1}-dimensional matrix containing the residuals.}
\item{sig}{Estimated standar deviation for the model.}
\item{fpca}{A list object containing model details, such as number of basis functions, number of principal components, and grid points used for the functional predictor variable.}
}
\author{Ufuk Beyaztas, Abhijit Mandal and Han Lin Shang}
\examples{
# Generate data
sim.dat <- data_generation(n = 100, j = 101, rho = 0.5, sig.e = 1, out.p = 0.05)
y <- sim.dat$y
x <- sim.dat$x
wei_mat <- sim.dat$w

flmPCA <- flm_pca(y=y, x=x)
}
