\name{bootstrap}
\alias{bootstrap}
\alias{print.resample.statistic}
\title{Generate Bootstrap Data and Statistics}
\description{
  Generates bootstrap samples from the original data and computes the
  bootstrap estimate of standard error and bias for a given statistic
  (if this is given).
}
\usage{
bootstrap(x, nb = 1, statistic = NULL, b = NULL,
          type = c("stationary","block"), \dots)
}
\arguments{
  \item{x}{a numeric vector or time series giving the original data.}
  \item{nb}{the number of bootstrap series to compute.}
  \item{statistic}{a function which when applied to a time series
    returns a vector containing the statistic(s) of interest.}
  \item{b}{if \code{type} is \code{"stationary"}, then \code{b} is the
    mean block length. If \code{type} is \code{"block"}, then \code{b}
    is the fixed block length.}
  \item{type}{the type of bootstrap to generate the simulated time
    series.  The possible input values are \code{"stationary"}
    (stationary bootstrap with mean block length \code{b}) and
    \code{"block"} (moving blocks bootstrap with block length
    \code{b}).}
  \item{\dots}{additional arguments for \code{statistic} which are
    passed unchanged each time \code{statistic} is called.}
}
\details{
  If \code{type} is \code{"stationary"}, then the stationary
  bootstrap scheme with mean block length \code{b} generates the
  simulated series. If \code{type} is \code{"block"}, then the moving
  blocks bootstrap with block length \code{b} generates the
  simulated series. 

  For consistency, the (mean) block length \code{b} should grow with
  \code{n} as \code{const * n^(1/3)}, where \code{n} is the number of
  observations in \code{x}. Note, that in general \code{const} depends
  on intricate properties of the process \code{x}. The default value for
  \code{const} has been determined by a Monte Carlo simulation using a
  Gaussian AR(1) (AR(1)-parameter of 0.5, 500 observations) process for
  \code{x}. It is chosen such that the mean square error for
  the bootstrap estimate of the variance of the empirical mean is
  minimized.

  Missing values are not allowed.

  There is a special print method for objects of class
  \code{"resample.statistic"} which by default uses
  \code{max(3, getOption("digits") - 3)} digits to format real numbers.
}
\value{
  If \code{statistic} is \code{NULL}, then it returns a matrix or time
  series with \code{nb} columns and \code{length(x)} rows containing the
  bootstrap data. Each column contains one bootstrap sample.

  If \code{statistic} is given, then a list of class
  \code{"resample.statistic"} with the following elements is returned:
  \item{statistic}{the results of applying \code{statistic} to each of
    the simulated time series.}
  \item{orig.statistic}{the results of applying \code{statistic} to the
    original series.}
  \item{bias}{the bias of the statistics computed as in a bootstrap
    setup.}
  \item{se}{the standard error of the statistics computed as in a
    bootstrap setup.}
  \item{call}{the original call of \code{bootstrap}.}
}
\references{
  H. R. Kuensch (1989):
  The Jackknife and the Bootstrap for General Stationary Observations.
  \emph{The Annals of Statistics} \bold{17}, 1217--1241.

  D. N. Politis and J. P. Romano (1994):
  The Stationary Bootstrap.
  \emph{Journal of the American Statistical Association} \bold{89},
  1303--1313.
}
\author{A. Trapletti}
\seealso{
  \code{\link{sample}},
  \code{\link{surrogate}}
}
\examples{
n <- 500  # Generate AR(1) process
e <- rnorm(n)  
x <- double(n)
x[1] <- rnorm(1)
for(i in 2:n) {
  x[i] <- 0.5 * x[i-1] + e[i]
}
x <- ts(x)

theta <- function(x)  # Autocorrelations up to lag 10
  return(acf(x, plot=FALSE)$acf[2:11])

bootstrap(x, nb=50, statistic=theta) 
}
\keyword{ts}
