\name{arma}
\title{Fit ARMA Models to Time Series}
\usage{
arma (x, order = c(1, 1), lag = NULL, coef = NULL, include.intercept = TRUE,
      series = NULL, qr.tol = 1e-07, ...)
coef (object)
residuals (object)
fitted (object)
print (object, digits = max(3,.Options$digits-3))
summary (object)
plot (object, ask = interactive())
print.summary (object, digits = max(3,.Options$digits-3),
               signif.stars = .Options$show.signif.stars, \dots)
}
\alias{arma}
\alias{coef.arma}
\alias{residuals.arma}
\alias{fitted.arma}
\alias{print.arma}
\alias{summary.arma}
\alias{plot.arma}
\alias{print.summary.arma}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{order}{a two dimensional integer vector giving the orders of the
    model to fit. \code{order[1]} corresponds to the AR part and
    \code{order[2]} to the MA part.}
  \item{lag}{a list with components \code{ar} and \code{ma}. Each
    component is an integer vector, specifying the AR and MA lags that are
    included in the model. If both, \code{order} and \code{lag}, are
    given, only the specification from \code{lag} is used.}
  \item{coef}{If given this numeric vector is used as the initial estimate
    of the ARMA coefficients. The preliminary estimator suggested in
    Hannan and Rissanen (1982) is used for the default initialization.}
  \item{include.intercept}{Should the model contain an intercept?}
  \item{series}{name for the series. Defaults to
    \code{deparse(substitute(x))}.}
  \item{qr.tol}{the \code{tol} argument for \code{\link{qr}} when computing
    the asymptotic standard errors of \code{coef}.}
  \item{\dots}{additional arguments for \code{\link{optim}} when fitting
    the model.}
  \item{object}{a fit from \code{arma}.}
  \item{digits, signif.stars}{see \code{\link{print.coefmat}}.}
  \item{ask}{Should the \code{plot} method work interactively? See
    \code{\link{interactive}}.}
  \item{\dots}{additional arguments for \code{\link{print}}.} 
}
\description{
  Fit an ARMA model to a univariate time series by conditional least
  squares. For exact maximum likelihood estimation see
  \code{\link{arima0}}. 
}
\details{
  The following parametrization is used for the ARMA(p,q) model:
  
  \code{y[t] = a[0] + a[1]y[t-1] + \dots + a[p]y[t-p] + b[1]e[t-1] +
    \dots + b[q]e[t-q] + e[t]},
  
  where \code{a[0]} is set to zero if no intercept is included. By using
  the argument \code{lag}, it is possible to fit a parsimonious submodel
  by setting arbitrary \code{a[i]} and \code{b[i]} to zero.
  
  \code{arma} uses \code{\link{optim}} to minimize the conditional
  sum-of-squared errors. The gradient is computed, if it is needed, by
  a finite-difference approximation. Default initialization is done by
  fitting a pure high-order AR model (see \code{\link{ar.ols}}). 
  The estimated residuals are then used for computing a least squares
  estimator of the full ARMA model. See Hannan and Rissanen (1982) for
  details.
  
  \code{summary} computes the asymptotic standard errors of the
  coefficient estimates from the numerically differentiated Hessian
  matrix approximation. The AIC is computed from the conditional
  sum-of-squared errors and not from the true maximum likelihood
  function. That may be problematic.
}
\value{
  For \code{arma} and its methods \code{print} and \code{plot} a list
  of class \code{"arma"} with the following elements:
  \item{lag}{the lag specification of the fitted model.}
  \item{coef}{estimated ARMA coefficients for the fitted model.}
  \item{css}{the conditional sum-of-squared errors.}
  \item{n.used}{the number of observations of \code{x}.}
  \item{residuals}{the series of residuals.}
  \item{fitted.values}{the fitted series.}
  \item{series}{the name of the series \code{x}.}
  \item{frequency}{the frequency of the series \code{x}.}
  \item{call}{the call of the \code{arma} function.}
  \item{asy.se.coef}{the asymptotic-theory standard errors of the
    coefficient estimates.}
  \item{convergence}{The \code{convergence} integer code from
    \code{\link{optim}}.}
  \item{include.intercept}{Does the model contain an intercept?}
  
  For \code{coef}, a numeric vector, for \code{residuals} and
  \code{fitted} a univariate time series.
  
  For \code{summary} and \code{print.summary} a list of
  class \code{"summary.arma"}.
}
\references{
  E. J. Hannan and J. Rissanen (1982): Recursive Estimation of Mixed
  Autoregressive-Moving Average Order. \emph{Biometrika} \bold{69},
  81-94. 
}
\author{
  A. Trapletti
}
\seealso{
  \code{\link{arima0}}, \code{\link{ar}}
}
\examples{
data (tcm)  
r <- diff(tcm10y)
summary (r.arma <- arma (r, order = c(1, 0)))
summary (r.arma <- arma (r, order = c(2, 0)))
summary (r.arma <- arma (r, order = c(0, 1)))
summary (r.arma <- arma (r, order = c(0, 2)))
summary (r.arma <- arma (r, order = c(1, 1)))
plot (r.arma)

data (nino)
s <- nino3.4
summary (s.arma <- arma (s, order=c(20,0)))
summary (s.arma
         <- arma (s, lag=list(ar=c(1,3,7,10,12,13,16,17,19),ma=NULL)))
acf (residuals(s.arma), na.action=na.remove)
pacf (residuals(s.arma), na.action=na.remove)
summary (s.arma
         <- arma (s, lag=list(ar=c(1,3,7,10,12,13,16,17,19),ma=12)))
summary (s.arma
         <- arma (s, lag=list(ar=c(1,3,7,10,12,13,16,17),ma=12)))
plot (s.arma)
}
\keyword{ts}


