\name{cross spectra}
\title{Cross Spectra Estimation}
\usage{
cross.spectrum (x, y, k = fejer.kernel(nextn(length(x))\%/\%2-1,length(x)/10))
co.spectrum (x, y, k = fejer.kernel(nextn(length(x))\%/\%2-1,length(x)/10), pl = TRUE, ...)
quadrature.spectrum (x, y, k = fejer.kernel(nextn(length(x))\%/\%2-1,length(x)/10), pl = TRUE, ...)
absolute.coherency (x, y, k = fejer.kernel(nextn(length(x))\%/\%2-1,length(x)/10), pl = TRUE, ...)
phase.spectrum (x, y, k = fejer.kernel(nextn(length(x))\%/\%2-1,length(x)/10), pl = TRUE, ...)
amplitude.spectrum (x, y, k = fejer.kernel(nextn(length(x))\%/\%2-1,length(x)/10), pl = TRUE, ...)
}
\alias{cross.spectrum}
\alias{co.spectrum}
\alias{quadrature.spectrum}
\alias{absolute.coherency}
\alias{phase.spectrum}
\alias{amplitude.spectrum}
\arguments{
  \item{x,y}{numeric vectors or time series.}
  \item{k}{a smoothing kernel of class \code{"kernel"}.}
  \item{pl}{a logical indicating whether the estimated spectrum is plotted.}
  \item{\dots}{additional arguments to \code{\link{plot.tsparam}}.}
}
\description{
  Estimates the different types of cross spectra of \code{x} and
  \code{y} by smoothing the bivariate periodogram with the kernel
  \code{k} (see Brockwell and Davies (1991), pp. 434-453). If \code{pl}
  is \code{TRUE}, then the estimated spectra
  are plotted. Additionally, the plot contains information about the
  "equivalent degrees of freedom" \code{df} and the bandwidth \code{bw}
  of \code{k} (equation (6.10) of T. Thrall (1983)). If \code{k} is a
  Daniell(0) kernel, then the periodogram itself is returned. For
  \code{absolute.coherency}, \code{phase.spectrum}, and
  \code{amplitude.spectrum} also the 95\%confidence bounds are plotted
  (see Brockwell and Davies (1991), equations (11.7.8), (11.7.9), and
  (11.7.11). Note, the error in (11.7.9)). 
  
  Missing values are not handled.
}
\value{
  A \code{"tsparam"} object.
}
\references{
  P. J. Brockwell and R. A. Davis (1991):
  \emph{Time Series: Theory and Methods}, 
  2nd Edition, Springer Verlag, NY, pp. 434-453.

  T. Thrall (1983): Computer Programming of spectrum estimation, in
  D. R. Brillinger and P. R. Krishnaiah (Eds.), \emph{Handbook of
  Statistics 3, Time Series in the Frequency Domain}, Elsevier Science,
  North-Holland, pp. 409-437.
}
\author{A. Trapletti}
\seealso{
  \code{\link{tsparam}}, \code{\link{kernel}}, \code{\link{spectrum}}
}
\examples{
data (sales)  # Reproduce example 11.7.1, Brockwell and Davis (1991)
sal <- diff (sales)
le <- diff (lead)
k <- daniell.kernel (6)
spectrum (le, k)
spectrum (sal, k)
absolute.coherency (le, sal, k, ylim = range(0,1))
phase.spectrum (le, sal, k)
}
\keyword{ts}

