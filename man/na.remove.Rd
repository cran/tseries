\name{na.remove}
\title{NA Handling Routines for Time Series}
\usage{
na.remove (x)
}
\alias{na.remove}
\alias{na.remove.ts}
\alias{na.remove.default}
\arguments{
  \item{x}{a numeric matrix, vector, univariate, or multivariate time series.} 
}
\description{
  Observations with missing values in some of the variables are
  removed.
}
\details{
  For \code{na.remove.ts} this changes the ``intrinsic'' time
  scale. It is assumed that both, the new and the old time scale are
  synchronized at the first and the last valid observation. In between,
  the new series is equally spaced in the new time scale. 
}
\value{
  An object without missing values. The attribute \code{"na.removed"}
  contains the indices of the removed missing values in \code{x}.
}
\author{A. Trapletti}
\seealso{\code{\link{na.omit}}, \code{\link{na.fail}}}
\examples{
x<-ts(c(5453.08,5409.24,5315.57,5270.53, # one and a half week stock index
        5211.66,NA,NA,5160.80,5172.37))  # data including a weekend
na.remove(x)  # eliminate weekend and introduce ``business'' time scale
}
\keyword{ts}
