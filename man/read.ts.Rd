\name{read.ts}
\alias{read.ts}
\title{Read Time Series Data}
\description{
  Reads a time series file.
}
\usage{
read.ts(file, header = FALSE, sep = "", skip = 0, \dots)
}
\arguments{
  \item{file}{the name of the file which the data are to be read from.
    Each line of the file contains one observation of the variables.}
  \item{header}{a logical value indicating whether the file contains
    the names of the variables as its first line.}
  \item{sep}{the field separator character.  Values on each line of the
    file are separated by this character.}
  \item{skip}{the number of lines of the data file to skip before
    beginning to read data.}
  \item{\dots}{Additional arguments for \code{\link{ts}} such as, e.g.,
    \code{start}.}
}
\details{
  Each row of the file represents an observation and each column
  contains a variable.  The first row possibly contains the names of the
  variables.
}
\author{A. Trapletti}
\seealso{
  \code{\link{ts}}.
}
\examples{
data(sunspots)
st <- start(sunspots)
fr <- frequency(sunspots)
write(sunspots, "sunspots", ncolumns=1)
x <- read.ts("sunspots", start=st, frequency=fr)
plot(x)
unlink("sunspots")
}
\keyword{file}
\keyword{ts}
