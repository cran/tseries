\name{acf}
\title{Autocorrelation Function}
\usage{
acf (x, lag = length(x)-1, correlation = TRUE, pl = TRUE, \dots)
}
\alias{acf}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{lag}{a scalar lag parameter.}
  \item{correlation}{a logical indicating whether the autocorrelation or
    autocovariance function is computed.}
  \item{pl}{a logical indicating whether the
    autocorrelation (covariance) function is plotted.}
  \item{\dots}{additional arguments to \code{\link{plot.tsparam}}.}
}
\description{
  Computes the sample autocorrelation (covariance) function of \code{x} up
  to lag \code{lag}. If \code{pl} is \code{TRUE}, then
  the autocorrelation (covariance) function is plotted. For the
  autocorrelation function also the 95\% confidence
  bounds for strict white noise are computed and plotted. Uses \code{fft}
  for efficiency reasons. 
  
  Missing values are not handled.
}
\value{
  A \code{"tsparam"} object.
}
\references{
  P. J. Brockwell and R. A. Davis (1991):
  \emph{Time Series: Theory and Methods}, 
  2nd Edition, Springer Verlag, NY, pp. 220-225, 373-375.
}
\author{A. Trapletti}
\seealso{
  \code{\link{tsparam}}
}
\examples{
data (sales)
acf (sales, lag = 20)
sal <- diff (sales)
acf (sal, lag = 20)
x <- rnorm (sal)
cr <- acf (x, lag = 20)
plot (cr, type = "o")
}
\keyword{ts}

