\name{ccf}
\title{Crosscorrelation Function}
\usage{
ccf (x, y, lag = length(x)-1, correlation = TRUE, pl = TRUE, \dots)
}
\alias{ccf}
\arguments{
  \item{x,y}{two numeric vector or time series.}
  \item{lag}{a scalar lag parameter.}
  item{correlation}{a logical indicating whether the crosscorrelation or
    crosscovariance function is computed.}
  \item{pl}{a logical indicating whether the
    crosscorrelation (covariance) function is plotted.}
  \item{\dots}{additional arguments to \code{\link{plot.tsparam}}.}
}
\description{
  Computes the sample crosscorrelation (covariance) function of \code{x} and
  \code{y} up to lag \code{lag}. If \code{pl} is \code{TRUE}, then
  the crosscorrelation (covariance) function is plotted. For the
  crosscorrelation function also the 95\% confidence
  bounds for strict white noise are plotted. Uses \code{fft}
  for efficiency reasons. 

  Missing values are not handled.
}
\value{
  A \code{"tsparam"} object.
}
\references{
  P. J. Brockwell and R. A. Davis (1991):
  \emph{Time Series: Theory and Methods}, 
  2nd Edition, Springer Verlag, NY, pp. 373-375, 405-417.
}
\author{A. Trapletti}
\seealso{
  \code{\link{tsparam}}
}
\examples{
data (sales)  # parts of Example 11.2.2 from Brockwell and Davies (1991).
sal <- diff (sales)
led <- diff(lead)
ccf (led, sal, lag = 20, ylim = range(-1,1), type="o")
}
\keyword{ts}



