\name{runs.test}
\alias{runs.test}
\title{Runs Test}
\description{
  Computes the runs test for randomness of the dichotomous (binary) data
  series \code{x}.
}
\usage{
runs.test(x, alternative = c("two.sided", "less", "greater"))
}
\arguments{
  \item{x}{a dichotomous factor.}
  \item{alternative}{indicates the alternative hypothesis and must be
    one of \code{"two.sided"} (default), \code{"less"}, or
    \code{"greater"}. You can specify just the initial letter.}
}
\details{
  This test searches for randomness in the observed data series
  \code{x} by examining the frequency of runs. A "run" is defined as a
  series of similar responses.

  Note, that by using the alternative \code{"less"} the null of
  randomness is tested against some kind of "under-mixing"
  ("trend"). By using the alternative \code{"greater"} the null of 
  randomness is tested against some kind of "over-mixing"
  ("mean-reversion").  
  
  Missing values are not allowed.
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{statistic}{the value of the test statistic.}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string indicating what type of test was
    performed.} 
  \item{data.name}{a character string giving the name of the data.}
  \item{alternative}{a character string describing the alternative
    hypothesis.}
}
\references{
  S. Siegel (1956): \emph{Nonparametric Statistics for the Behavioural
  Sciences}, McGraw-Hill, New York.
  
  S. Siegel and N. J. Castellan (1988): \emph{Nonparametric Statistics for the
  Behavioural Sciences}, 2nd edn, McGraw-Hill, New York.
}
\author{A. Trapletti}
\examples{
x <- factor(sign(rnorm(100)))  # randomness
runs.test(x)

x <- factor(rep(c(-1,1),50))  # over-mixing
runs.test(x)
}
\keyword{ts}
