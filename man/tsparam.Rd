\name{tsparam}
\title{Time Series Paramater Object}
\usage{
tsparam (x, y, xlab="", ylab="", lab="", ybnd=NULL, ysgnf=NULL)
plot (obj, xlim = range(obj$x), ylim = range(attr(obj,"ybnd"),obj$y),
      xlab = attr(obj,"xlab"), ylab = attr(obj,"ylab"),
      type = "l", main = "", grid = TRUE, col = "black", 
      bcol = "blue", gridcol = "lightgrey") 
print (obj, digits = max(3,.Options$digits-3), \dots)
}
\alias{tsparam}
\alias{plot.tsparam}
\alias{print.tsparam}
\arguments{
  \item{x,y}{the vectors containg the data.}
  \item{xlab,ylab}{the labels of \code{x,y}.}
  \item{lab}{an additional label which contains a more detailed description of the
    \code{"tsparam"} object or the way it was build.}
  \item{ybnd}{a list of vectors with the same length as \code{x} and
    \code{y} which contains confidence bounds for the \code{y} data.}
  \item{ysgnf}{a vector of characters with the same length as \code{x}
    containing the significance codes for the \code{y} data.}
  \item{obj}{a \code{"tsparam"} object.}
  \item{xlim,ylim}{the ranges of the axes.}
  \item{xlab,ylab}{the axis labels.}
  \item{type}{the type of the plot.}
  \item{main}{a main title string.}
  \item{grid}{a logical indicating whether a grid is plotted.}
  \item{col, bcol, gridcol}{strings specifying the colors for
    the data, the confidence bounds, and the grid.}
  \item{digits}{the number of digits to format real numbers.}
  \item{\dots}{additional arguments to \code{\link{print}}.}
}
\description{
  The \code{"tsparam"} class is designed to represent results from time
  series analysis. It just provides common \code{\link{plot}} and 
  \code{\link{print}} \code{\link{methods}} for estimated parameters.
}
\value{
  \code{tsparam} constructs and returns a class \code{"tsparam"} object.
  It consists of a list with elements \code{x} and \code{y}. Attached
  attributes are \code{"xlab"}, \code{"ylab"}, \code{"ybnd"},
  \code{"ysgnf"}, and \code{"lab"}.
  
  The \code{print} and \code{plot} methods return the invisible input
  object.
}
\author{A. Trapletti}
\seealso{
  \code{\link{plot}}, \code{\link{print}}
}
\examples{
x <- 1:100
y <- rnorm (100)
t <- tsparam (x,y,"index","random numbers")
plot (t)
print (t)
}
\keyword{ts}

