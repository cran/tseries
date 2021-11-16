\name{bds.test}
\alias{bds.test}
\alias{print.bdstest}
\title{BDS Test}
\description{
  Computes and prints the BDS test statistic for the null that \code{x}
  is a series of i.i.d. random variables.
}
\usage{
bds.test(x, m = 3, eps = seq(0.5 * sd(x), 2 * sd(x), length.out = 4),
         trace = FALSE)
}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{m}{an integer indicating that the BDS test statistic is computed
    for embedding dimensions \code{2}, \dots, \code{m}.}
  \item{eps}{a numeric vector of epsilon values for close points.  The
    BDS test is computed for each element of \code{eps}.  It should be
    set in terms of the standard deviation of \code{x}.}
  \item{trace}{a logical indicating whether some informational output is
    traced.}
}
\details{
  This test examines the ``spatial dependence'' of the observed
  series.  To do this, the series is embedded in \code{m}-space and the
  dependence of \code{x} is examined by counting ``near'' points.
  Points for which the distance is less than \code{eps} are called
  ``near''.  The BDS test statistic is asymptotically standard Normal.
  
  Missing values are not allowed.

  There is a special print method for objects of class \code{"bdstest"}
  which by default uses 4 digits to format real numbers.
}
\value{
  A list with class \code{"bdstest"} containing the following components:
  \item{statistic}{the values of the test statistic.}
  \item{p.value}{the p-values of the test.}
  \item{method}{a character string indicating what type of test was
    performed.}
  \item{parameter}{a list with the components \code{m} and \code{eps}
    containing the embedding dimensions and epsilon values for which the
    statistic is computed.}
  \item{data.name}{a character string giving the name of the data.}
}
\references{
  J. B. Cromwell, W. C. Labys and M. Terraza (1994):
  \emph{Univariate Tests for Time Series Models},
  Sage, Thousand Oaks, CA, pages 32--36.
}
\author{B. LeBaron, Ported to R by A. Trapletti}
\examples{
x <- rnorm(100)
bds.test(x)  # i.i.d. example

x <- c(rnorm(50), runif(50))
bds.test(x)  # not identically distributed

x <- quadmap(xi = 0.2, a = 4.0, n = 100)
bds.test(x)  # not independent
}
\keyword{ts}
