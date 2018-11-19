\name{kpss.test}
\alias{kpss.test}
\title{KPSS Test for Stationarity}
\description{
  Computes the Kwiatkowski-Phillips-Schmidt-Shin (KPSS) test for the
  null hypothesis that \code{x} is level or trend stationary.
}
\usage{
kpss.test(x, null = c("Level", "Trend"), lshort = TRUE)
}
\arguments{
  \item{x}{a numeric vector or univariate time series.}
  \item{null}{indicates the null hypothesis and must be one of
    \code{"Level"} (default) or \code{"Trend"}.  You can specify just
    the initial letter.}
  \item{lshort}{a logical indicating whether the short or long version
    of the truncation lag parameter is used.}
}
\details{
  To estimate \code{sigma^2} the Newey-West estimator is used.
  If \code{lshort} is \code{TRUE}, then the truncation lag parameter is
  set to \code{trunc(4*(n/100)^0.25)}, otherwise
  \code{trunc(12*(n/100)^0.25)} is used.  The p-values are interpolated
  from Table 1 of Kwiatkowski et al. (1992). If the computed statistic
  is outside the table of critical values, then a warning message is
  generated.
  
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
}
\references{
  D. Kwiatkowski, P. C. B. Phillips, P. Schmidt, and Y. Shin (1992):
  Testing the Null Hypothesis of Stationarity against the Alternative of
  a Unit Root.
  \emph{Journal of Econometrics} \bold{54}, 159--178.
}
\author{A. Trapletti}
\seealso{
  \code{\link{pp.test}}
}
\examples{
x <- rnorm(1000)  # is level stationary
kpss.test(x)

y <- cumsum(x)  # has unit root
kpss.test(y)

x <- 0.3*(1:1000)+rnorm(1000)  # is trend stationary
kpss.test(x, null = "Trend")
}
\keyword{ts}
