\name{adf.test}
\alias{adf.test}
\title{Augmented Dickey--Fuller Test}
\description{
  Computes the Augmented Dickey-Fuller test for the null that \code{x} has
  a unit root.
}
\usage{
adf.test(x, alternative = c("stationary", "explosive"),
         k = trunc((length(x)-1)^(1/3)))
}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{alternative}{indicates the alternative hypothesis and must be
    one of \code{"stationary"} (default) or \code{"explosive"}. You can
    specify just the initial letter.}
  \item{k}{the lag order to calculate the test statistic.}
}
\details{
  The general regression equation which incorporates a constant and a
  linear trend is used and the t-statistic for a first order
  autoregressive coefficient equals one is computed. The number of lags
  used in the regression is \code{k}. The default value of
  \code{trunc((length(x)-1)^(1/3))} corresponds to the suggested upper
  bound on the rate at which the number of lags, \code{k}, should be
  made to grow with the sample size for the general \code{ARMA(p,q)}
  setup. Note that for \code{k} equals zero the standard Dickey-Fuller
  test is computed. The p-values are interpolated from Table 4.2, p. 103
  of Banerjee et al. (1993). If the computed statistic is outside the
  table of critical values, then a warning message is generated. 
  
  Missing values are not allowed.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{the lag order.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name of the data.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
}
\references{
  A. Banerjee, J. J. Dolado, J. W. Galbraith, and D. F. Hendry (1993):
  \emph{Cointegration, Error Correction, and the Econometric Analysis of
    Non-Stationary Data},
  Oxford University Press, Oxford.
  
  S. E. Said and D. A. Dickey (1984):
  Testing for Unit Roots in Autoregressive-Moving Average Models of
  Unknown Order.
  \emph{Biometrika} \bold{71}, 599--607.
}
\author{A. Trapletti}
\seealso{
  \code{\link{pp.test}}
}
\examples{
x <- rnorm(1000)  # no unit-root
adf.test(x)

y <- diffinv(x)   # contains a unit-root
adf.test(y)
}
\keyword{ts}
