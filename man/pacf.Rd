\name{pacf}
\title{Partial Autocorrelation Function}
\usage{
pacf (x, lag = length(x) - 1, pl = TRUE, \dots)
}
\alias{pacf}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{lag}{a scalar lag parameter.}
  \item{pl}{a logical indicating whether the
    partial autocorrelation function is plotted.}
  \item{\dots}{additional arguments to \code{\link{plot.tsparam}}.}
}
\description{
  Computes the sample partial autocorrelation function of \code{x} up
  to lag \code{lag}. If \code{pl} is \code{TRUE}, then
  the partial autocorrelation function and the 95\% confidence
  bounds for strict white noise are also plotted. 
  
  Missing values are not handled.
}
\details{
  The partial autocorrelations are obtained from Yule-Walker estimates
  of the successive autoregressive processes. The algorithm of Durbin
  is used (see Box and Jenkins (1976), p. 65).
}
\value{
  A \code{"tsparam"} object.
}
\references{
  P. J. Brockwell and R. A. Davis (1991):
  \emph{Time Series: Theory and Methods}, 
  2nd Edition, Springer Verlag, NY, pp. 98-102.

  G. E. P. Box and G. M. Jenkins (1976):
  \emph{Time Series Analysis, Forecasting and Control},
  Holden-Day, San Francisco, pp. 64, 65.
}
\author{A. Trapletti}
\seealso{
  \code{\link{tsparam}}
}
\examples{
n <- 1500
e <- rnorm (n, sd=2.25^2)  # Example 8.2.1 from Brockwell and Davies (1991)
x <- double (n)
x[1] <- 0
x[2] <- 0
for (i in 3:n)
{
  x[i] <- 0.5*x[i-1]+0.2*x[i-2]+e[i]
}
x <- window(ts(x),start=501)
plot (x)
acf (x, 40, type="o", ylim=range(-1,1))
pacf (x, 40, type="o", ylim=range(-1,1))
}
\keyword{ts}

