\name{print.bdstest}
\title{Print BDS Test}
\usage{
print (x, digits = 4)
}
\alias{print.bdstest}
\arguments{
  \item{x}{a list with class \code{"bdstest"}.}
  \item{digits}{the number of digits to format real numbers.}
}
\description{
  Generic function for printing lists with class \code{"bdstest"}.
}
\author{A. Trapletti}
\seealso{
  \code{\link{print}}, \code{\link{bds.test}}.
}
\examples{
x <- rnorm (100)
obj <- bds.test(x)
print (obj)
}
\keyword{ts}

