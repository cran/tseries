\name{surrogate}
\alias{surrogate}
\title{Generate Surrogate Data and Statistics}
\description{
  Generates \code{ns} surrogate samples from the original data \code{x}
  and computes the standard error and the bias of \code{statistic} as in
  a bootstrap setup, if \code{statistic} is given.
}
\usage{
surrogate(x, ns = 1, fft = FALSE, amplitude = FALSE,
          statistic = NULL, \dots)
}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{ns}{the number of surrogate series to compute.}
  \item{fft}{a logical indicating whether phase randomized surrogate
    data is generated.}
  \item{amplitude}{a logical indicating whether amplitude-adjusted
    surrogate data is computed.}
  \item{statistic}{a function which when applied to a time series
    returns a vector containing the statistic(s) of interest.}
  \item{\dots}{Additional arguments for \code{statistic} which are
    passed unchanged each time it is called.}
}
\details{
  If \code{fft} is \code{FALSE}, then \code{x} is mixed in temporal
  order, so that all temporal dependencies are eliminated, but the
  histogram of the original data is preserved.  If \code{fft} is
  \code{TRUE}, then surrogate data with the same spectrum as \code{x} is
  computed by randomizing the phases of the Fourier coefficients of
  \code{x}.  If in addition \code{amplitude} is \code{TRUE}, then also
  the amplitude distribution of the original series is preserved.

  Note, that the interpretation of the computed standard error and bias
  is different than in a bootstrap setup.

  To compute the phase randomized surrogate and the amplitude adjusted
  data algorithm 1 and 2 from Theiler et al. (1992), pp. 183, 184 are
  used.
  
  Missing values are not allowed.
}
\value{
  If \code{statistic} is \code{NULL}, then it returns a matrix or time
  series with \code{ns} columns and \code{length(x)} rows containing the
  surrogate data. Each column contains one surrogate sample.

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
  \item{call}{the original call of \code{surrogate}.}
}
\references{
  J. Theiler, B. Galdrikian, A. Longtin, S. Eubank, and J. D. Farmer
  (1992): Using Surrogate Data to Detect
  Nonlinearity in Time Series, in \emph{Nonlinear Modelling and
  Forecasting}, Eds. M. Casdagli and S. Eubank, Santa Fe Institute,
  Addison Wesley, 163--188.
}
\author{A. Trapletti}
\seealso{
  \code{\link{sample}}, \code{\link{tsbootstrap}}
}
\examples{
x <- 1:10  # Simple example
surrogate(x)

n <- 500  # Generate AR(1) process
e <- rnorm(n)  
x <- double(n)
x[1] <- rnorm(1)
for(i in 2:n) {
  x[i] <- 0.4 * x[i-1] + e[i]
}
x <- ts(x)

theta <- function(x)  # Autocorrelations up to lag 10
  return(acf(x, plot=FALSE)$acf[2:11])

surrogate(x, ns=50, fft=TRUE, statistic=theta) 
}
\keyword{ts}
