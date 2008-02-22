\name{arma}
\alias{arma}
\title{Fit ARMA Models to Time Series}
\description{
  Fit an ARMA model to a univariate time series by conditional least
  squares.  For exact maximum likelihood estimation see
  \code{\link{arima0}}. 
}
\usage{
arma(x, order = c(1, 1), lag = NULL, coef = NULL,
     include.intercept = TRUE, series = NULL, qr.tol = 1e-07, \dots)
}
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
}
\details{
  The following parametrization is used for the ARMA(p,q) model:
  
  \deqn{y[t] = a[0] + a[1]y[t-1] + \dots + a[p]y[t-p] + b[1]e[t-1] +
    \dots + b[q]e[t-q] + e[t],}
  
  where \eqn{a[0]} is set to zero if no intercept is included. By using
  the argument \code{lag}, it is possible to fit a parsimonious submodel
  by setting arbitrary \eqn{a[i]} and \eqn{b[i]} to zero.
  
  \code{arma} uses \code{\link{optim}} to minimize the conditional
  sum-of-squared errors. The gradient is computed, if it is needed, by
  a finite-difference approximation. Default initialization is done by
  fitting a pure high-order AR model (see \code{\link{ar.ols}}). 
  The estimated residuals are then used for computing a least squares
  estimator of the full ARMA model. See Hannan and Rissanen (1982) for
  details.
}
\value{
  A list of class \code{"arma"} with the following elements:
  \item{lag}{the lag specification of the fitted model.}
  \item{coef}{estimated ARMA coefficients for the fitted model.}
  \item{css}{the conditional sum-of-squared errors.}
  \item{n.used}{the number of observations of \code{x}.}
  \item{residuals}{the series of residuals.}
  \item{fitted.values}{the fitted series.}
  \item{series}{the name of the series \code{x}.}
  \item{frequency}{the frequency of the series \code{x}.}
  \item{call}{the call of the \code{arma} function.}
  \item{vcov}{estimate of the asymptotic-theory covariance matrix for the
    coefficient estimates.}
  \item{convergence}{The \code{convergence} integer code from
    \code{\link{optim}}.}
  \item{include.intercept}{Does the model contain an intercept?}
}
\references{
  E. J. Hannan and J. Rissanen (1982):
  Recursive Estimation of Mixed Autoregressive-Moving Average
  Order.
  \emph{Biometrika} \bold{69}, 81--94. 
}
\author{
  A. Trapletti
}
\seealso{
  \code{\link{summary.arma}} for summarizing ARMA model fits;
  \code{\link{arma-methods}} for further methods;
  \code{\link{arima0}}, \code{\link{ar}}.
}
\examples{
data(tcm)  
r <- diff(tcm10y)
summary(r.arma <- arma(r, order = c(1, 0)))
summary(r.arma <- arma(r, order = c(2, 0)))
summary(r.arma <- arma(r, order = c(0, 1)))
summary(r.arma <- arma(r, order = c(0, 2)))
summary(r.arma <- arma(r, order = c(1, 1)))
plot(r.arma)

data(nino)
s <- nino3.4
summary(s.arma <- arma(s, order=c(20,0)))
summary(s.arma
         <- arma(s, lag=list(ar=c(1,3,7,10,12,13,16,17,19),ma=NULL)))
acf(residuals(s.arma), na.action=na.remove)
pacf(residuals(s.arma), na.action=na.remove)
summary(s.arma
         <- arma(s, lag=list(ar=c(1,3,7,10,12,13,16,17,19),ma=12)))
summary(s.arma
         <- arma(s, lag=list(ar=c(1,3,7,10,12,13,16,17),ma=12)))
plot(s.arma)
}
\keyword{ts}
