\name{cumulative.periodogram}
\title{Cumulative Periodogram}
\usage{
cumulative.periodogram (x, pl = TRUE, ...)
}
\alias{cumulative.periodogram}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{pl}{a logical indicating whether the cumulative periodogram is plotted.}
  \item{\dots}{additional arguments to \code{\link{plot.tsparam}}.}
}
\description{
  Computes the cumulative periodogram and the Kolmogorov-Smirnov test
  for hidden periodicities of the series \code{x}. (see Brockwell and
  Davies (1991), pp. 334-342). If \code{pl} is \code{TRUE}, then the
  cumulative periodogram together with the 95\%confidence bounds for
  Gaussian white noise are plotted.
    
  Missing values are not handled.
}
\value{
  A \code{"tsparam"} object.
}
\references{
  P. J. Brockwell and R. A. Davis (1991):
  \emph{Time Series: Theory and Methods}, 
  2nd Edition, Springer Verlag, NY, pp. 334-342. 
}
\author{A. Trapletti}
\seealso{
  \code{\link{tsparam}}, \code{\link{spectrum}}
}
\examples{
data (sales)  
sal <- diff (sales)
cumulative.periodogram (sal)
x <- rnorm(sal)
cumulative.periodogram (x)
}
\keyword{ts}

