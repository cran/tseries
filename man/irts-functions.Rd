\name{irts-functions}
\alias{irts-functions}
\alias{approx.irts}
\alias{daysecond}
\alias{is.businessday}
\alias{is.weekend}
\alias{read.irts}
\alias{weekday}
\alias{write.irts}
\title{Basic Functions for Irregular Time-Series Objects}
\description{
  Basic functions related to irregular time-series objects.
}
\usage{
daysecond(object, tz = "GMT")
approx.irts(object, time, \dots)
is.businessday(object, tz = "GMT")
is.weekend(object, tz = "GMT")
read.irts(file, format = "\%Y-\%m-\%d \%H:\%M:\%S", tz = "GMT", \dots)
weekday(object, tz = "GMT")
write.irts(object, file = "", append = FALSE, quote = FALSE,
           sep = " ", eol = "\n", na = "NA", dec = ".",
           row.names = FALSE, col.names = FALSE, qmethod = "escape",
           format = "\%Y-\%m-\%d \%H:\%M:\%S", tz = "GMT", usetz = FALSE,
           format.value = NULL, \dots)
}
\arguments{
  \item{object}{an object of class \code{"irts"}; usually, a result 
    of a call to \code{\link{irts}}.}
  \item{format, tz, usetz}{formatting related arguments, see
    \code{\link{format.POSIXct}}.} 
  \item{time}{an object of class \code{"POSIXct"} specifying the times
    at which to interpolate the irregularly spaced time-series.}
  \item{file, append, quote, sep, eol, na, dec, row.names, col.names,
    qmethod}{reading and writing related arguments, see
    \code{\link{read.table}} and \code{\link{write.table}}.}
  \item{format.value}{a string which specifies the formatting of the
    values when writing an irregular time-series object to a
    file. \code{format.value} is passed unchanged as argument
    \code{format} to the function \code{\link{formatC}}.} 
  \item{\dots}{further arguments passed to or from other methods: for
    \code{approx.irts} passed to \code{\link{approx}}; for
    \code{read.irts} passed to \code{\link{read.table}}; for
    \code{write.irts} passed to \code{\link{data.frame}}.}
}
\details{
  \code{daysecond} and \code{weekday} return the number of seconds since
  midnight (the same day) and the weekday as a decimal number (0-6,
  Sunday is 0), respectively.

  \code{is.businessday} and \code{is.weekend} test which entries of an
  irregular time-series object are recorded on business days and
  weekends, respectively. 

  \code{approx.irts} interpolates an irregularly spaced time-series at
  prespecified times.

  \code{read.irts} is the function to read irregular time-series
  objects from a file.
  
  \code{write.irts} is the function to write irregular time-series
  objects to a file.
}
\value{
  For \code{daysecond} and \code{weekday} a vector of decimal numbers
  representing the number of seconds and the weekday, respectively.

  For \code{is.businessday} and \code{is.weekend} a vector of
  \code{"logical"} representing the test results for each time.

  For \code{approx.irts}, \code{read.irts} and \code{write.irts} an
  object of class \code{"irts"}.
}
\author{
  A. Trapletti
}
\seealso{
  \code{\link{irts}},
  \code{\link{irts-methods}}
}
\examples{
n <- 10
t <- cumsum(rexp(n, rate = 0.1))
v <- rnorm(n)
x <- irts(t, v)

daysecond(x)
weekday(x)
is.businessday(x)
is.weekend(x)
x

approx.irts(x, seq(ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "GMT"),
            by = "10 secs", length.out = 7), rule = 2)

\dontrun{
file <- tempfile()

# To write an irregular time-series object to a file one might use
write.irts(x, file = file)

# To read an irregular time-series object from a file one might use
read.irts(file = file)

unlink(file)
}
}
\keyword{ts}
