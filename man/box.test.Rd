\name{box.test}
\title{Box-Pierce and Box-Ljung Test}
\usage{
box.test (x, lag = 1, pierce = TRUE) 
}
\alias{box.test}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{lag}{a scalar lag parameter.}
  \item{pierce}{a logical indicating whether the
    Box-Pierce or the Box-Ljung test is computed.}
}
\description{
  If \code{pierce} is \code{TRUE}, then the Box-Pierce test
  for examining the null of independence in the time series \code{x} is
  computed. Else the Box-Ljung statistic is computed. Uses
  \code{lag} autocorrelation coefficients for the statistic.
  
  Missing values are not handled.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{the degrees of freedom of the approximate chi-square
    distribution of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name of the data.}
}
\references{
  A. C. Harvey (1993): \emph{Time Series Models},
  2nd Edition, Harvester Wheatsheaf, NY, pp. 44, 45.
}
\author{A. Trapletti}
\examples{
x <- rnorm (100)
box.test (x, lag = 1, pierce = TRUE)
}
\keyword{ts}

