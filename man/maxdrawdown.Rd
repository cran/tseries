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
  A list containing the following three components:
  \item{maxdrawdown}{double representing the max drawdown or max
    loss statistic.}
  \item{from}{the index (or vector of indices) where the max drawdown
    period starts.} 
  \item{to}{the index (or vector of indices) where the max drawdown
    period ends.} 
}
\author{A. Trapletti}
\seealso{
  \code{\link{sterling}}
}
\examples{
# Toy example
x <- c(1:10, 9:7, 8:14, 13:8, 9:20)
mdd <- maxdrawdown(x)
mdd

plot(x)
segments(mdd$from, x[mdd$from], mdd$to, x[mdd$from], col="grey")
segments(mdd$from, x[mdd$to], mdd$to, x[mdd$to], col="grey")
mid <- (mdd$from + mdd$to)/2
arrows(mid, x[mdd$from], mid, x[mdd$to], col="red", length = 0.16)

# Realistic example
data(EuStockMarkets)
dax <- log(EuStockMarkets[,"DAX"])
mdd <- maxdrawdown(dax)
mdd

plot(dax)
segments(time(dax)[mdd$from], dax[mdd$from],
         time(dax)[mdd$to], dax[mdd$from], col="grey")
segments(time(dax)[mdd$from], dax[mdd$to],
         time(dax)[mdd$to], dax[mdd$to], col="grey")
mid <- time(dax)[(mdd$from + mdd$to)/2]
arrows(mid, dax[mdd$from], mid, dax[mdd$to], col="red", length = 0.16)
}
\keyword{ts}
