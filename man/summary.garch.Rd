\name{summary.garch}
\alias{summary.garch}
\alias{print.summary.garch}
\title{Summarizing GARCH Model Fits}
\description{
  Methods for creating and printing summaries of GARCH model fits.
}
\usage{
\method{summary}{garch}(object, \dots)
\method{print}{summary.garch}(x, digits = max(3, getOption("digits") - 3),
      signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"garch"}; usually, a result of a
    call to \code{\link{garch}}.}
  \item{x}{an object of class \code{"summary.garch"}; usually, a result
    of a call to the summary method for objects of class \code{"garch"}.}
  \item{digits, signif.stars}{see \code{\link{printCoefmat}}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  \code{summary} computes the asymptotic standard errors of the
  coefficient estimates from an outer-product approximation of the
  Hessian evaluated at the estimates, see Bollerslev (1986).  It
  furthermore tests the residuals for normality and remaining ARCH
  effects, see \code{\link{jarque.bera.test}} and
  \code{\link{Box.test}}. 
}
\value{
  A list of class \code{"summary.garch"}.
}
\references{
  T. Bollerslev (1986):
  Generalized Autoregressive Conditional Heteroscedasticity.
  \emph{Journal of Econometrics} \bold{31}, 307--327.
}
\seealso{
  \code{\link{garch}}
}
\keyword{models}
\keyword{ts}
