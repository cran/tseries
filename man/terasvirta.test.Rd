\name{terasvirta.test}
\title{Teraesvirta Neural Network Test for Nonlinearity}
\usage{
terasvirta.test (x, lag = 1, type = c("chisq","F"), scale = TRUE)
terasvirta.test (x, y, type = c("chisq","F"), scale = TRUE)
}
\alias{terasvirta.test}
\alias{terasvirta.test.ts}
\alias{terasvirta.test.default}
\arguments{
  \item{x}{a numeric vector, matrix, or time series.}
  \item{y}{a numeric vector.}
  \item{lag}{an integer which specifies the model order in terms of
    lags.} 
  \item{type}{a string indicating whether the Chi-Square test or the
    F-test is computed. Valid types are \code{"chisq"} and \code{"F"}.}
  \item{scale}{a logical indicating whether the data should be scaled 
    before computing the test statistic. The default arguments to
    \code{\link{scale}} are used.}
}
\description{
  Generically computes Teraesvirta's neural network test for neglected
  nonlinearity either for the time series \code{x} or the regression
  \code{y~x}.
}
\details{
  The null is the hypotheses of linearity in
  ``mean''. This test uses a Taylor series expansion of the activation
  function to arrive at a suitable test statistic. If \code{type} equals
  \code{"F"}, then the F-statistic instead of the Chi-Square statistic
  is used in analogy to the classical linear regression. 
  
  Missing values are not allowed.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{parameter}{a list containing the additional parameters used to
    compute the test statistic.}
  \item{data.name}{a character string giving the name of the data.}
  \item{arguments}{additional arguments used to compute the test statistic.}
}
\references{
  T. Teraesvirta, C. F. Lin, and C. W. J. Granger (1993): Power of the
  Neural Network Linearity Test. \emph{Journal of Time Series Analysis}
  14, 209-220.
}
\author{A. Trapletti}
\seealso{
  \code{\link{white.test}}
}
\examples{
n <- 1000

x <- runif (1000, -1, 1)  # Non-linear in ``mean'' regression 
y <- x^2 - x^3 + 0.1*rnorm(x)
terasvirta.test (x, y)

terasvirta.test (cbind(x,x^2,x^3), y)  # Is the polynomial of order 2 misspecified?

x[1] <- 0.0
for (i in (2:n))  # Generate time series which is nonlinear in ``mean''
{
  x[i] <- 0.4*x[i-1] + tanh(x[i-1]) + rnorm (1, sd=0.5)
}
x <- as.ts(x)
plot (x)
terasvirta.test (x)
}
\keyword{ts}




