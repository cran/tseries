\name{read.matrix}
\alias{read.matrix}
\title{Read Matrix Data}
\description{
  Reads a matrix data file.
}
\usage{
read.matrix(file, header = FALSE, sep = "", skip = 0)
}
\arguments{
  \item{file}{the name of the file which the data are to be read from.}
  \item{header}{a logical value indicating whether the file contains the
    names of the columns as its first line.}
  \item{sep}{the field separator character.  Values on each line of the
    file are separated by this character.}
  \item{skip}{the number of lines of the data file to skip before
    beginning to read data.}
}
\details{
  Usually each row of the file represents an observation and each column
  contains a variable.  The first row possibly contains the names of the
  variables (columns).

  \code{read.matrix} might be more efficient than
  \code{\link{read.table}} for very large data sets.
}
\author{A. Trapletti}
\seealso{
  \code{\link{read.table}}.
}
\examples{
x <- matrix(0, 10, 10)
write(x, "test", ncolumns=10)
x <- read.matrix("test")
x
unlink("test")
}
\keyword{file}
\keyword{ts}
