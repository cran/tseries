\name{sharpe}
\alias{sharpe}
\title{Sharpe Ratio}
\description{
  This function computes the Sharpe ratio of the univariate time series
  (or vector) \code{x}.  
}
\usage{
sharpe(x, r = 0, scale = sqrt(250))
}
\arguments{
  \item{x}{a numeric vector or univariate time series corresponding to a
    portfolio's cumulated returns.}
  \item{r}{the risk free rate. Default corresponds to using portfolio
    returns not in excess of the riskless return.} 
  \item{scale}{a scale factor. Default corresponds to an annualization
    when working with daily financial time series data.}  
}
\details{
  The Sharpe ratio is defined as a portfolio's mean return in excess of
  the riskless return divided by the portfolio's standard deviation. In
  finance the Sharpe Ratio represents a measure of the portfolio's
  risk-adjusted (excess) return. 
}
\value{
  a double representing the Sharpe ratio.
}
\author{A. Trapletti}
\seealso{
  \code{\link{sterling}}
}
\examples{
data(EuStockMarkets)
dax <- log(EuStockMarkets[,"DAX"])
ftse <- log(EuStockMarkets[,"FTSE"])
sharpe(dax)
sharpe(ftse)
}
\keyword{ts}
