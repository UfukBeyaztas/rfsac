\name{data_generation}
\alias{data_generation}
\title{Generate a dataset for the spatial scalar-on-function regression model}
\description{This function can be used to generate a dataset for the spatial scalar-on-function regression model \deqn{
Y = \rho W Y + \int X(t) \beta(t) dt + \epsilon,} where \eqn{Y} denotes the scalar response, \eqn{X(t)} denotes the functional predictor, \eqn{W} is the spatial weight matrix, \eqn{\beta(t)} denotes the regression coefficient function, \eqn{\rho} is the spatial autocorrelation parameter, and \eqn{\epsilon} is the error process.}
\usage{
data_generation(n, j, rho, sig.e, out.p = 0)
}
\arguments{
  \item{n}{An integer, specifying the number of observations for each variable to be generated.}
  \item{j}{An integer, denoting the number of grid points, i.e., a fine grid on the interval [0, 1].}
  \item{rho}{A numeric value, denoting the spatial autocorrelation parameter.}
  \item{sig.e}{A numeric value, denoting the standard deviation of the error process.}
  \item{out.p}{A numeric value between 0 and 1 specifying the outlier percentage.}
}
\details{In the data generation process, first, the  functional predictor is generated based on the following process: \deqn{
X(t) = \sum_{j=1}^5 \kappa_j v_j(t),} where \eqn{ \kappa_j } is a vector generated from a Normal distribution with mean zero and variance \eqn{4 j^{-3/2}} and \deqn{v_j(t) = \sin(j \pi t) - \cos(j \pi t).} The  regression coefficient function is generated from \eqn{\sin(2 \pi t)}. The error process is generated from the standard normal distribution. The row-normalized spatial weight matrix is constructed based on a one-dimensional regular grid, with \eqn{w_{i i'} = \frac{1}{d_{i i'}}} where \eqn{d_{i i'} = \vert i - i' \vert}, for \eqn{i \neq i'}.}
\value{A list object with the following components:
\item{y}{An \eqn{n \times 1}-dimensional matrix containing the observations of simulated scalar response variable.}
\item{x}{A matrix with dimension \eqn{n \times j} containing the observations of simulated functional predictor variable.}
\item{w}{A matrix with dimension \eqn{n \times n}, i.e., the generated spatial weight matrix.}
\item{tcoefs}{A vector with length \eqn{n} containing the generated regression coefficient function.}
\item{out.index}{A vector consisting of the index set of outlying observations.}
}
\author{Ufuk Beyaztas, Abhijit Mandal and Han Lin Shang}
\examples{
# Generate data
sim.dat <- data_generation(n = 100, j = 101, rho = 0.5, sig.e = 1, out.p = 0.05)

y <- sim.dat$y
x <- sim.dat$x
wei_mat <- sim.dat$w

# true coefficient
btrue <- c(sim.dat$tcoefs$data)

library(fda.usc)
fx <- fdata(x, argvals = seq(0, 1, length.out = 101))
par(mfrow = c(1,2))
plot(y, type = "p", pch = 16, xlab = "Index", ylab = "", main = "Response")
plot(fx, lty = 1, ylab = "", xlab = "Grid point",
main = expression(X(t)), mgp = c(2, 0.5, 0))
dev.off()
}

