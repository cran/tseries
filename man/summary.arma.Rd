\name{summary.arma}
\alias{summary.arma}
\alias{print.summary.arma}
\title{Summarizing ARMA Model Fits}
\description{
  Methods for creating and printing summaries of ARMA model fits.
}
\usage{
\method{summary}{arma}(object, \dots)
\method{print}{summary.arma}(x, digits = max(3, getOption("digits") - 3),
      signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"arma"}; usually, a result of a
    call to \code{\link{arma}}.}
  \item{x}{an object of class \code{"summary.arma"}; usually, a result
    of a call to the summary method for objects of class \code{"arma"}.}
  \item{digits, signif.stars}{see \code{\link{printCoefmat}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  The summary method computes the asymptotic standard errors of the
  coefficient estimates from the numerically differentiated Hessian
  matrix approximation.  The AIC is computed from the conditional
  sum-of-squared errors and not from the true maximum likelihood
  function. That may be problematic.
}
\value{
  A list of class \code{"summary.arma"}.
}
\seealso{
  \code{\link{arma}}
}
\keyword{models}
\keyword{ts}
