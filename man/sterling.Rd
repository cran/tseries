\name{sterling}
\alias{sterling}
\title{Sterling Ratio}
\description{
  This function computes the Sterling ratio of the univariate time series
  (or vector) \code{x}.  
}
\usage{
sterling(x)
}
\arguments{
  \item{x}{a numeric vector or univariate time series corresponding to a
    portfolio's cumulated returns.}
}
\details{
  The Sterling ratio is defined as a portfolio's overall return divided
  by the portfolio's \code{\link{maxdrawdown}} statistic. In finance the
  Sterling Ratio represents a measure of the portfolio's risk-adjusted
  return.   
}
\value{
  a double representing the Sterling ratio.
}
\author{A. Trapletti}
\seealso{
  \code{\link{maxdrawdown}}, \code{\link{sharpe}}
}
\examples{
data(EuStockMarkets)
dax <- log(EuStockMarkets[,"DAX"])
ftse <- log(EuStockMarkets[,"FTSE"])
sterling(dax)
sterling(ftse)
}
\keyword{ts}
