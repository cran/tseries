\name{spectrum}
\title{Spectrum Estimation}
\usage{
spectrum (x, k = fejer.kernel(nextn(length(x))\%/\%2-1,length(x)/10), pl = TRUE, ...)
}
\alias{spectrum}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{k}{a smoothing kernel of class \code{"kernel"}.}
  \item{pl}{a logical indicating whether the estimated spectrum is plotted.}
  \item{\dots}{additional arguments to \code{\link{plot.tsparam}}.}
}
\description{
  Computes the "discrete spectral average estimator" of the spectrum of
  \code{x}, i.e., the periodogram of \code{x} is smoothed with the given weights
  from the discrete kernel \code{k} (see Brockwell and Davies (1991), p. 351, for a
  definition). If \code{pl} is \code{TRUE}, then the estimated spectrum
  is plotted. Additionally, the plot contains information about the
  "equivalent degrees of freedom" \code{df} and the bandwidth \code{bw}
  of \code{k} (equation (6.10) of T. Thrall (1983)). If \code{k} is a
  Daniell(0) kernel, then the periodogram
  itself is returned.
  
  Missing values are not handled.
}
\value{
  A \code{"tsparam"} object.
}
\references{
  P. J. Brockwell and R. A. Davis (1991):
  \emph{Time Series: Theory and Methods}, 
  2nd Edition, Springer Verlag, NY, pp. 350-365.

  T. Thrall (1983): Computer Programming of spectrum estimation, in
  D. R. Brillinger and P. R. Krishnaiah (Eds.), \emph{Handbook of
  Statistics 3, Time Series in the Frequency Domain}, Elsevier Science,
  North-Holland, pp. 409-437.
}
\author{A. Trapletti}
\seealso{
  \code{\link{tsparam}}, \code{\link{kernel}}
}
\examples{
data (sales)  # Reproduce parts of example 11.7.1, Brockwell and Davis (1991)
sal <- diff (sales)
le <- diff (lead)
k <- daniell.kernel (6)
spectrum (le, k)
spectrum (sal, k)
}
\keyword{ts}

