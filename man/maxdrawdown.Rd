\name{maxdrawdown}
\alias{maxdrawdown}
\title{Maximum Drawdown or Maximum Loss}
\description{
  This function computes the maximum drawdown or maximum loss of the
  univariate time series (or vector) \code{x}. 
}
\usage{
maxdrawdown(x)
}
\arguments{
  \item{x}{a numeric vector or univariate time series.}
}
\details{
  The max drawdown or max loss statistic is defined as the maximum
  value drop after one of the peaks of \code{x}. For financial
  instruments the max drawdown represents the worst investment loss for
  a buy-and-hold strategy invested in \code{x}.
}
\value{
  A double representing the max drawdown or max loss statistic.
}
\author{A. Trapletti}
\examples{
x <- c(1:10, 9:7, 8:14, 13:8, 9:20)
maxdrawdown(x)

plot(x)
lines(x = c(20, 28), y = c(14, 14), col="grey")
lines(x = c(20, 28), y = c(8, 8), col="grey")
arrows(26, 14, 26, 8, , col="red", length = 0.16)
}
\keyword{ts}
