\name{irts}
\alias{irts}
\alias{is.irts}
\alias{as.irts}
\alias{as.irts.default}
\alias{as.irts.zoo}
\title{Irregularly Spaced Time-Series}
\description{
  The function \code{irts} is used to create irregular time-series
  objects.  

  \code{as.irts} coerces an object to an irregularly spaced
  time-series. \code{is.irts} tests whether an object is an irregularly
  spaced time series. 
}
\usage{
irts(time, value)
as.irts(object)
is.irts(object)
}
\arguments{
  \item{time}{a numeric vector or a vector of class \code{"POSIXct"}
    representing the time-stamps of the irregular time-series
    object. The elements of the numeric vector are construed as the
    (signed) number of seconds since the beginning of 1970, see
    \code{\link{POSIXct}}.} 
  \item{value}{a numeric vector or matrix representing the values of the
    irregular time-series object.}
  \item{object}{an R object to be coerced to an irregular time-series
    object or an R object to be tested whether it is an irregular
    time-series object.} 
}
\details{
  The function \code{irts} is used to create irregular time-series
  objects. These are scalar or vector valued time series indexed by a
  time-stamp of class \code{"POSIXct"}. Unlike objects of class
  \code{"ts"}, they can be used to represent irregularly spaced
  time-series.
}
\value{
  A list of class \code{"irts"} with the following elements:
  \item{time}{a vector of class \code{"POSIXct"}.}
  \item{value}{a numeric vector or matrix.}
}
\author{
  A. Trapletti
}
\seealso{
  \code{\link{ts}},
  \code{\link{POSIXct}},
  \code{\link{irts-methods}},
  \code{\link{irts-functions}}
}
\examples{
n <- 10
t <- cumsum(rexp(n, rate = 0.1))
v <- rnorm(n)
x <- irts(t, v)
x

as.irts(cbind(t, v))

is.irts(x)

# Multivariate
u <- rnorm(n)
irts(t, cbind(u, v))
}
\keyword{ts}
