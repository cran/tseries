\name{intgrt}
\title{Discrete Integrals}
\usage{
intgrt <- function (x, lag = 1, differences = 1,
                    xi = rep(0.0,lag*differences*dim(as.matrix(x))[2]))
}
\alias{intgrt}
\alias{intgrt.ts}
\alias{intgrt.default}
\arguments{
  \item{x}{a numeric vector, matrix, or time series.}
  \item{lag}{a scalar lag parameter.}
  \item{differences}{an integer representing the order of the difference.}
  \item{xi}{a numeric vector, matrix, or time series containing the
    initial values for the integrals.}
}
\description{
  Computes the inverse function of the lagged differences
  \code{\link{diff}}. 

  Missing values are not handled.
}
\details{
  \code{intgrt} is a generic function with methods for \code{"ts"}
  and \code{default} for vector and matrix.
}
\value{
  A numeric vector, matrix, or time series representing the discrete integral of
  \code{x}. 
}
\author{A. Trapletti}
\seealso{
  \code{\link{diff}}
}
\examples{
s <- 1:10
d <- diff(s)
intgrt(d,xi=1)
}
\keyword{ts}

