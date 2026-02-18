\name{jarque.bera.test}
\alias{jarque.bera.test}
\title{Jarque--Bera Test}
\description{
  Tests the null of normality for \code{x} using the Jarque-Bera
  test statistic.
}
\usage{
jarque.bera.test(x)
}
\arguments{
  \item{x}{a numeric vector or time series.}
}
\details{
  This test is a joint statistic using skewness and kurtosis
  coefficients.
  
  Missing values are not allowed.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{parameter}{the degrees of freedom.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name of the data.}
}
\references{
  \bibinfo{Cromwell+Labys+Terraza:1994}{note}{Pages 20--22}
  \bibshow{Cromwell+Labys+Terraza:1994}
}
\author{A. Trapletti}
\examples{
x <- rnorm(100)  # null
jarque.bera.test(x)

x <- runif(100)  # alternative
jarque.bera.test(x)
}
\keyword{ts}
