\name{plotOHLC}
\alias{plotOHLC}
\title{Plot Open--High--Low--Close Bar Chart}
\description{
  Plots open--high--low--close bar chart of a (financial) time series.
}
\usage{
plotOHLC(x, xlim = NULL, ylim = NULL, xlab = "Time", ylab, col = par("col"),
         bg = par("bg"), axes = TRUE, frame.plot = axes, ann = par("ann"),
         main = NULL, date = c("calendar", "julian"), format = "\%Y-\%m-\%d", \dots)
}
\arguments{
  \item{x}{a multivariate time series object of class \code{"mts"}.}
  \item{xlim, ylim}
  \item{xlab, ylab}
  \item{col, bg}
  \item{axes, frame.plot}
  \item{ann, main}{graphical arguments, see \code{\link{plot}},
    \code{\link{plot.default}} and \code{\link{par}}.}
  \item{date}{a string indicating the type of x axis
    annotation. Default is calendar dates.}
  \item{format}{a string indicating the format of the x axis annotation if
    \code{date == "calendar"}. For details see
    \code{\link{format.POSIXct}}.} 
  \item{\dots}{further graphical arguments passed to
    \code{\link{plot.window}}, \code{\link{title}}, \code{\link{axis}},
    and \code{\link{box}}.}
}
\details{Within an open--high--low--close bar chart, each bar represents
  price information for the time interval between the open and the close
  price. The left tick for each bar indicates the open price for the
  time interval. The right tick indicates the closing price for the time
  interval. The vertical length of the bar represents the price range
  for the time interval.
}
\author{A. Trapletti}
\seealso{
  \code{\link{plot.default}},
  \code{\link{format.POSIXct}},
  \code{\link{get.hist.quote}}
}
\examples{
\dontrun{
# Plot OHLC bar chart for the last `nDays' days of the instrument `instrument'

nDays <- 50
instrument <- "^spc"

start <- strftime(as.POSIXlt(Sys.time() - nDays*24*3600), format="\%Y-\%m-\%d") 
end <- strftime(as.POSIXlt(Sys.time()), format = "\%Y-\%m-\%d") 
x <- get.hist.quote(instrument = instrument, start = start, end = end)

plotOHLC(x, ylab = "price", main = instrument)
}
}
\keyword{hplot}
\keyword{ts}
