\name{white.test}
\title{White Neural Network Test for Nonlinearity}
\alias{white.test}
\alias{white.test.ts}
\alias{white.test.default}
\description{
  Generically computes the White neural network test for neglected
  nonlinearity either for the time series \code{x} or the regression
  \code{y~x}.
}
\usage{
\method{white.test}{ts}(x, lag = 1, qstar = 2, q = 10, range = 4, 
            type = c("Chisq","F"), scale = TRUE, \dots)
\method{white.test}{default}(x, y, qstar = 2, q = 10, range = 4,
            type = c("Chisq","F"), scale = TRUE, \dots)
}
\arguments{
  \item{x}{a numeric vector, matrix, or time series.}
  \item{y}{a numeric vector.}
  \item{lag}{an integer which specifies the model order in terms of
    lags.} 
  \item{q}{an integer representing the number of phantom hidden units
    used to compute the test statistic.}
  \item{qstar}{the test is conducted using \code{qstar} principal
    components of the phantom hidden units. The first principal
    component is omitted since in most cases it appears to be collinear
    with the input vector of lagged variables. This strategy preserves
    power while still conserving degrees of freedom.}
  \item{range}{the input to hidden unit weights are initialized uniformly
    over [-range/2, range/2].}
  \item{type}{a string indicating whether the Chi-Squared test or the
    F-test is computed. Valid types are \code{"Chisq"} and \code{"F"}.}
  \item{scale}{a logical indicating whether the data should be scaled 
    before computing the test statistic. The default arguments to
    \code{\link{scale}} are used.}
  \item{\dots}{further arguments to be passed from or to methods.}
}
\details{
  The null is the hypotheses of linearity in ``mean''. This
  type of test is consistent against arbitrary nonlinearity
  in mean. If \code{type} equals \code{"F"}, then the F-statistic
  instead of the Chi-Squared statistic is used in analogy to the
  classical linear regression. 
  
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
  T. H. Lee, H. White, and C. W. J. Granger (1993): Testing for
  neglected nonlinearity in time series models. \emph{Journal of
  Econometrics} \bold{56}, 269-290.
}
\author{A. Trapletti}
\seealso{
  \code{\link{terasvirta.test}}
}
\examples{
n <- 1000

x <- runif(1000, -1, 1)  # Non-linear in ``mean'' regression 
y <- x^2 - x^3 + 0.1*rnorm(x)
white.test(x, y)

## Is the polynomial of order 2 misspecified?
white.test(cbind(x,x^2,x^3), y)

## Generate time series which is nonlinear in ``mean''
x[1] <- 0.0
for(i in (2:n)) {
  x[i] <- 0.4*x[i-1] + tanh(x[i-1]) + rnorm(1, sd=0.5)
}
x <- as.ts(x)
plot(x)
white.test(x)
}
\keyword{ts}
