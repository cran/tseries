\name{pp.test}
\alias{pp.test}
\title{Phillips--Perron Unit Root Test}
\description{
  Computes the Phillips-Perron test for the null hypothesis that
  \code{x} has a unit root.
}
\usage{
pp.test(x, alternative = c("stationary", "explosive"),
        type = c("Z(alpha)", "Z(t_alpha)"), lshort = TRUE)
}
\arguments{
  \item{x}{a numeric vector or univariate time series.}
  \item{alternative}{indicates the alternative hypothesis and must be
    one of \code{"stationary"} (default) or \code{"explosive"}.  You can
    specify just the initial letter.}
  \item{type}{indicates which variant of the test is computed and must
    be one of \code{"Z(alpha)"} (default) or \code{"Z(t_alpha)"}.}
  \item{lshort}{a logical indicating whether the short or long version
    of the truncation lag parameter is used.}
}
\details{
  The general regression equation which incorporates a constant and a
  linear trend is used and the \code{Z(alpha)} or \code{Z(t_alpha)}
  statistic for a first order autoregressive coefficient equals one are
  computed.  To estimate \code{sigma^2} the Newey-West estimator is
  used.  If \code{lshort} is \code{TRUE}, then the truncation lag
  parameter is set to \code{trunc(4*(n/100)^0.25)}, otherwise
  \code{trunc(12*(n/100)^0.25)} is used.  The p-values are interpolated
  from Table 4.1 and 4.2, p. 103 of Banerjee et al. (1993).  If the
  computed statistic is outside the table of critical values, then a
  warning message is generated.
  
  Missing values are not handled.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{the truncation lag parameter.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name of the data.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
}
\references{
  A. Banerjee, J. J. Dolado, J. W. Galbraith, and D. F. Hendry (1993):
  \emph{Cointegration, Error Correction, and the Econometric Analysis
    of Non-Stationary Data}, Oxford University Press, Oxford.
  
  P. Perron (1988):
  Trends and Random Walks in Macroeconomic Time Series.
  \emph{Journal of Economic Dynamics and Control} \bold{12}, 297--332.
}
\author{A. Trapletti}
\seealso{
  \code{\link{adf.test}}
}
\examples{
x <- rnorm(1000)  # no unit-root
pp.test(x)

y <- cumsum(x)  # has unit root
pp.test(y)
}
\keyword{ts}
