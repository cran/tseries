\name{adf.test}
\title{Augmented Dickey-Fuller Test}
\usage{
adf.test (x, k = trunc((length(x)-1)^(1/3)))
}
\alias{adf.test}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{k}{the lag order to calculate the test statistic.}
}
\description{
  Computes the Augmented Dickey-Fuller test for the null that \code{x} has
  a unit root. The general regression
  equation which incorporates a constant and a linear trend is used and
  the t-statistic for a first order autoregressive coefficient
  equals one is computed. The number of lags used in the regression is
  \code{k}. The default value of \code{trunc((length(x)-1)^(1/3))}
  corresponds to the suggested upper bound on the rate at which the
  number of lags, \code{k}, should be made to grow with the sample size
  for the general \code{ARMA(p,q)} setup. Note that for \code{k} equals
  zero the standard Dickey-Fuller test is computed. The p-values are
  interpolated from Table 4.2, p. 103 of Banerjee et al. (1993).
  
  Missing values are not handled.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{the lag order.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name of the data.}
}
\references{
  A. Banerjee, J. J. Dolado, J. W. Galbraith, and D. F. Hendry (1993):
  \emph{Cointegration, Error Correction, and the Econometric Analysis of
  Non-Stationary Data}, Oxford University Press, Oxford.
  
  S. E. Said and D. A. Dickey (1984): Testing for Unit Roots in
  Autoregressive-Moving Average Models of Unknown
  Order. \emph{Biometrika} 71, 599-607.
}
\author{A. Trapletti}
\examples{
x <- rnorm (1000)
adf.test (x)
y <- intgrt (x)
adf.test (y)
}
\keyword{ts}

