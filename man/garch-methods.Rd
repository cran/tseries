\name{garch-methods}
\alias{garch-methods}
\alias{predict.garch}
\alias{coef.garch}
\alias{residuals.garch}
\alias{fitted.garch}
\alias{print.garch}
\alias{plot.garch}
\alias{logLik.garch}
\title{Methods for Fitted GARCH Models}
\description{
  Methods for fitted GARCH model objects, typically the result of a call
  to \code{garch}.
}
\usage{
\method{predict}{garch}(object, newdata, genuine = FALSE, \dots)
\method{coef}{garch}(object, \dots)
\method{residuals}{garch}(object, \dots)
\method{fitted}{garch}(object, \dots)
\method{print}{garch}(x, digits = max(3, getOption("digits") - 3), \dots)
\method{plot}{garch}(x, ask = interactive(), \dots)
\method{logLik}{garch}(object, \dots)
}
\arguments{
  \item{object, x}{an object of class \code{"garch"}; usually, a result
    of a call to \code{\link{garch}}.}
  \item{newdata}{a numeric vector or time series to compute GARCH
    predictons.  Defaults to \code{eval(parse(text=object$series))}.}
  \item{genuine}{a logical indicating whether a genuine prediction
    should be made, i.e., a prediction for which there is no target
    observation available.}
  \item{digits, signif.stars}{see \code{\link{print.coefmat}}.}
  \item{ask}{Should the \code{plot} method work interactively?  See
    \code{\link{interactive}}.}
  \item{\dots}{further arguments passed to or from other methods.} 
}
\details{
  \code{predict} returns +/- the conditional standard deviation
  predictions from a fitted GARCH model.

  \code{coef} returns the coefficient estimates.

  \code{residuals} returns the GARCH residuals, i.e., the time series
  used to fit the model divided by the computed conditional standard
  deviation predictions for this series. Under the assumption of
  conditional normality the residual series should be i.i.d. standard
  normal.  

  \code{fitted} returns +/- the conditional standard deviation
  predictions for the series which has been used to fit the model.

  \code{plot} graphically investigates normality and remaining ARCH
  effects for the residuals.

  \code{logLik} returns the log-likelihood value of the GARCH(p, q)
  model represented by \code{object} evaluated at the estimated
  coefficients. It is assumed that first max(p, q) values are fixed.
}  
\value{
  For \code{predict} a bivariate time series (two-column matrix) of
  predictions. 
  
  For \code{coef}, a numeric vector, for \code{residuals} and
  \code{fitted} a univariate (vector) and a bivariate time series
  (two-column matrix), respectively.

  For \code{plot} and \code{print}, the fitted GARCH model object.
}
\keyword{models}
\keyword{ts}
