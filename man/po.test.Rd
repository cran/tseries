\name{po.test}
\alias{po.test}
\title{Phillips--Ouliaris Cointegration Test}
\description{
  Computes the Phillips-Ouliaris test for the null hypothesis that
  \code{x} is not cointegrated.
}
\usage{
po.test(x, demean = TRUE, lshort = TRUE) 
}
\arguments{
  \item{x}{a matrix or multivariate time series.}
  \item{demean}{a logical indicating whether an intercept is included in
    the cointegration regression or not.}
  \item{lshort}{a logical indicating whether the short or long version
    of the truncation lag parameter is used.}
}
\details{
  The Phillips-Perron Z(alpha) statistic for a unit root in the
  residuals of the cointegration regression is computed, see also
  \code{\link{pp.test}}. The unit root is estimated from a regression of
  the first variable (column) of \code{x} on the remaining variables of
  \code{x} without a constant and a linear trend.  To estimate
  \code{sigma^2} the Newey-West estimator is used.  If \code{lshort} is
  \code{TRUE}, then the truncation lag parameter is set to
  \code{trunc(n/100)}, otherwise \code{trunc(n/30)} is used.  The
  p-values are interpolated from Table Ia and Ib, p. 189 of Phillips and
  Ouliaris (1990).  If the computed statistic is outside the table of
  critical values, then a warning message is generated.
  
  The dimension of \code{x} is restricted to six variables.  Missing
  values are not handled.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{the truncation lag parameter.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name of the data.}
}
\references{
  P. C. B. Phillips and S. Ouliaris (1990):
  Asymptotic Properties of Residual Based Tests for Cointegration.
  \emph{Econometrica} \bold{58}, 165--193.
}
\author{A. Trapletti}
\seealso{
  \code{\link{pp.test}}
}
\examples{
x <- ts(diffinv(matrix(rnorm(2000),1000,2)))  # no cointegration
po.test(x)

x <- diffinv(rnorm(1000))
y <- 2.0-3.0*x+rnorm(x,sd=5)
z <- ts(cbind(x,y))  # cointegrated
po.test(z)
}
\keyword{ts}
