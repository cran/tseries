\name{runs.test}
\alias{runs.test}
\title{Runs Test}
\description{
  Computes the runs test statistic for the null of independence for
  the series \code{x}.
}
\usage{
runs.test(x)
}
\arguments{
  \item{x}{a numeric vector or time series.}
}
\details{
  This test searches for independence in the observed series by
  examining the frequency of runs. Before computing the statistic,
  zero values of \code{x} are eliminated.
  
  Missing values are not allowed.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name of the data.}
}
\references{
  J. B. Cromwell, W. C. Labys and M. Terraza (1994):
  \emph{Univariate Tests for Time Series Models},
  Sage, Thousand Oaks, CA, 28--30.
}
\author{A. Trapletti}
\examples{
x <- rnorm(100)  # no runs
runs.test(x)

x <- 1:100  # one big run
runs.test(x)
}
\keyword{ts}
