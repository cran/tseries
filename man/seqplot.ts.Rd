\name{seqplot.ts}
\alias{seqplot.ts}
\title{Plot Two Time Series}
\description{
  Plot two time series on the same plot frame.
}
\usage{
seqplot.ts(x, y, colx = "black", coly = "red", typex = "l",
           typey = "l", pchx = 1, pchy = 1, ltyx = "solid",
           ltyy = "solid", oma = c(6, 0, 5, 0), ann = par("ann"),
           xlab = "Time", ylab = deparse(substitute(x)), main = NULL)
}

\arguments{
  \item{x, y}{the time series.}
  \item{colx, coly}{color code or name for the \code{x} and \code{y}
    series, see \code{\link{colors}}, \code{\link{palette}}.}
  \item{typex, typey}{what type of plot should be drawn for the \code{x}
    and \code{y} series, see \code{\link{plot}}.}
  \item{pchx, pchy}{character or integer code for kind of points/lines
    for the \code{x} and \code{y} series.}
  \item{ltyx, ltyy}{line type code for the \code{x} and \code{y} series,
    see \code{\link{lines}}.}
  \item{oma}{a vector giving the size of the outer margins in lines of
    text, see \code{\link{par}}.}
  \item{ann}{annotate the plots? See \code{\link{par}}.}
  \item{xlab, ylab}{titles for the x and y axis.}
  \item{main}{an overall title for the plot.}
}
\details{
  Unlike \code{\link{plot.ts}} the series can have different time bases,
  but they should have the same frequency. Unlike \code{\link{ts.plot}}
  the series can be plotted in different styles and for multivariate
  \code{x} and \code{y} the common variables are plotted together in a
  separate array element.
}
\value{
  None.
}
\author{A. Trapletti}
\seealso{
  \code{\link{ts}}, \code{\link{plot.ts}}
}
\examples{
data(USeconomic)
x <- ts.union(log(M1), log(GNP), rs, rl)
m.ar <- ar(x, method = "ols", order.max = 5)
y <- predict(m.ar, x, n.ahead = 200, se.fit = FALSE)
seqplot.ts(x, y)
}
\keyword{hplot}
\keyword{ts}
