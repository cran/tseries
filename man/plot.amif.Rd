\name{plot.amif}
\alias{plot.amif}
\title{Plot Method for Auto Mutual Information Functions}
\description{
  The \code{plot} method for auto mutual information function objects.
}
\usage{
\method{plot}{amif}(x, ci.col = "blue", ylab = "AMIF", \dots)
}
\arguments{
  \item{x}{an object of class \code{"amif"}, usually the result of a
    call to \code{\link{amif}}.}
  \item{ci.col}{color to plot the confidence bound.}
  \item{ylab}{the title for the y axis.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The AMIF is plotted along with a confidence bound for the null of
  surrogate data as computed by Monte-Carlo simulation percentiles.
}
\seealso{
  \code{\link{amif}}.
}
\keyword{hplot}
\keyword{ts}
