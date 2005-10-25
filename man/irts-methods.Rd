\name{irts-methods}
\alias{irts-methods}
\alias{lines.irts}
\alias{plot.irts}
\alias{points.irts}
\alias{print.irts}
\alias{time.irts}
\alias{value}
\alias{value.irts}
\alias{[.irts}
\title{Methods for Irregular Time-Series Objects}
\description{
  Methods for irregular time-series objects.
}
\usage{
\method{lines}{irts}(x, type = "l", \dots)
\method{plot}{irts}(x, type = "l", plot.type = c("multiple", "single"), 
     xlab = "Time", ylab = NULL, main = NULL, ylim = NULL,
     oma = c(6, 0, 5, 0), \dots)
\method{points}{irts}(x, type = "p", \dots)                    
\method{print}{irts}(x, format = "\%Y-\%m-\%d \%H:\%M:\%S", tz = "GMT",
      usetz = TRUE, format.value = NULL, \dots)
\method{time}{irts}(x, \dots)
\method{value}{irts}(x, \dots)
\method{[}{irts}(x, i, j, \dots)
}
\arguments{
  \item{x}{an object of class \code{"irts"}; usually, a result 
    of a call to \code{\link{irts}}.}
  \item{type, plot.type, xlab, ylab, main, ylim, oma}{graphical
    arguments, see \code{\link{plot}}, \code{\link{points}},
    \code{\link{lines}}, \code{\link{par}}, and \code{\link{plot.ts}}.}
  \item{format, tz, usetz}{formatting related arguments, see
    \code{\link{format.POSIXct}}.}
  \item{format.value}{a string which specifies the formatting of the
    values when printing an irregular time-series 
    object. \code{format.value} is passed unchanged as argument
    \code{format} to the function \code{\link{formatC}}.}
  \item{i, j}{indices specifying the parts to extract from an irregular
    time-series object.}
  \item{\dots}{further arguments passed to or from other methods: for
    \code{lines} passed to \code{\link{lines}}; for \code{plot} passed
    to \code{\link{plot}}, \code{\link{plot.default}}, and
    \code{\link{mtext}}; for \code{points} passed to
    \code{\link{points}}; for \code{print} passed to
    \code{\link{formatC}}; for \code{time}, \code{value}, and
    \code{[.irts} unused.} 
}
\details{
  \code{plot} is the method for plotting irregular time-series objects.

  \code{points} and \code{lines} are the methods for drawing a sequence
  of points as given by an irregular time-series object and joining the
  corresponding points with line segments, respectively.

  \code{print} is the method for printing irregular time-series objects.

  \code{time} and \code{value} are the methods for extracting the
  sequence of times and the sequence of values of an irregular
  time-series object.
  
  \code{[.irts} is the method for extracting parts of irregular
  time-series objects. 
}
\value{
  For \code{time} an object of class \code{"POSIXct"} representing the
  sequence of times. For \code{value} a vector or matrix representing
  the sequence of values. 

  For \code{[.irts} an object of class \code{"irts"} representing the
  extracted part.
  
  For \code{plot}, \code{points}, \code{lines}, and \code{print} the
  irregular time-series object. 
}
\author{
  A. Trapletti
}
\seealso{
  \code{\link{irts}},
  \code{\link{irts-functions}}
}
\examples{
n <- 10
t <- cumsum(rexp(n, rate = 0.1))
v <- rnorm(n)
x <- irts(t, v)

x
time(x)
value(x)
plot(x)
points(x)

t <- cumsum(c(t[1], rexp(n-1, rate = 0.2)))
v <- rnorm(n, sd = 0.1)
x <- irts(t, v)

lines(x, col = "red")
points(x, col = "red")

# Multivariate
t <- cumsum(rexp(n, rate = 0.1))
u <- rnorm(n)
v <- rnorm(n)
x <- irts(t, cbind(u, v))

x
x[,1]
x[1:3,]
x[1:3,1]
plot(x)
}
\keyword{ts}
