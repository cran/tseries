\name{surrogate}
\title{Generate Surrogate Data}
\usage{
surrogate (x, ns = 1, fft = FALSE, amplitude = FALSE)
}
\alias{surrogate}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{ns}{the number of surrogate series to compute.}
  \item{fft}{a logical indicating whether phase randomized surrogate
    data is generated.}
  \item{amplitude}{a logical indicating whether
    amplitude-adjusted surrogate data is computed.}
}
\description{
  Computes \code{ns} surrogate samples from the original data
  \code{x}. If \code{fft} is \code{FALSE}, then \code{x} is mixed in
  temporal order, so that all temporal dependencies are eliminated, but
  the histogram of the original data is preserved. If \code{fft}, then
  surrogate data with the same spectrum as \code{x} is computed by
  randomizing the phases of the Fourier coefficients of \code{x}. If in
  addition \code{amplitude} is \code{TRUE}, then also the amplitude
  distribution of the original series is preserved.
  
  Missing values are not allowed.
}
\value{
  A matrix with \code{ns} columns and \code{length{x}} rows containing
  the surrogate data. Each column contains one surrogate sample.
}
\references{
  J. Theiler, B. Galdrikian, A. Longtin, S. Eubank, and J. D. Farmer
  (1992): Using Surrogate Data to Detect
  Nonlinearity in Time Series, in \emph{Nonlinear Modelling and
  Forecasting}, Eds. M. Casdagli and S. Eubank, Santa Fe Institute,
  Addison Wesley, pp. 163-188.
}
\details{
  To compute the phase randomized surrogate and the amplitude adjusted
  data algorithm 1 and 2 from Theiler et al. (1992), pp. 183, 184 are
  used. 
}
\author{A. Trapletti}
\seealso{
  \code{\link{sample}}
}
\examples{
x <- 1:10
surrogate (x)
}
\keyword{ts}


