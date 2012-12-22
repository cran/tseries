\name{portfolio.optim}
\alias{portfolio.optim}
\alias{portfolio.optim.ts}
\alias{portfolio.optim.default}
\title{Portfolio Optimization}
\description{
  Computes an efficient portfolio from the given return series \code{x}
  in the mean-variance sense.
}
\usage{
\method{portfolio.optim}{default}(x, pm = mean(x), riskless = FALSE,
                shorts = FALSE, rf = 0.0, reslow = NULL, reshigh = NULL,
                covmat = cov(x), \dots)
}
\arguments{
  \item{x}{a numeric matrix or multivariate time series consisting of a
    series of returns.}
  \item{pm}{the desired mean portfolio return.}
  \item{riskless}{a logical indicating whether there is a riskless
    lending and borrowing rate.}
  \item{shorts}{a logical indicating whether shortsales on the risky
    securities are allowed.} 
  \item{rf}{the riskfree interest rate.}
  \item{reslow}{a vector specifying the (optional) lower bound on
    allowed portfolio weights.}  
  \item{reshigh}{a vector specifying the (optional) upper bound on
    allowed portfolio weights.} 
  \item{covmat}{the covariance matrix of asset returns.}
  \item{\dots}{further arguments to be passed from or to methods.}
}
\details{
  The computed portfolio has the desired expected return \code{pm} and
  no other portfolio exists, which has the same mean return, but a
  smaller variance. Inequality restrictions of the form \eqn{w_l \le w
  \le w_h}{w_l <= w <= w_h} can be imposed using the \code{reslow} and
  \code{reshigh} vectors. An alternative covariance matrix estimate can
  be supplied via the \code{covmat} argument. To solve the quadratic
  program, \code{\link[quadprog]{solve.QP}} is used.

  \code{portfolio.optim} is a generic function with methods for
  multivariate \code{"ts"} and \code{default} for matrix.

  Missing values are not allowed.
}
\value{
  A list containing the following components:
  \item{pw}{the portfolio weights.}
  \item{px}{the returns of the overall portfolio.}
  \item{pm}{the expected portfolio return.}
  \item{ps}{the standard deviation of the portfolio returns.}
}
\references{
  E. J. Elton and M. J. Gruber (1991): \emph{Modern Portfolio Theory and
  Investment Analysis}, 4th Edition, Wiley, NY, pp. 65-93.

  C. Huang and R. H. Litzenberger (1988): \emph{Foundations for
  Financial Economics}, Elsevier, NY, pp. 59-82.
}
\author{A. Trapletti}
\seealso{
  \code{\link[quadprog]{solve.QP}}
}
\examples{
x <- rnorm(1000)
dim(x) <- c(500,2)
res <- portfolio.optim(x)
res$pw

require("zoo")			# For diff() method.
X <- diff(log(as.zoo(EuStockMarkets)))
res <- portfolio.optim(X)                 ## Long only
res$pw
res <- portfolio.optim(X, shorts=TRUE)    ## Long/Short
res$pw
}
\keyword{ts}
