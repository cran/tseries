\name{garch}
\alias{garch}
\alias{garch.control}
\title{Fit GARCH Models to Time Series}
\description{
  Fit a Generalized Autoregressive Conditional Heteroscedastic GARCH(p,
  q) time series model to the data by computing the maximum-likelihood
  estimates of the conditionally normal model.
}
\usage{
garch(x, order = c(1, 1), series = NULL, control = garch.control(\dots), \dots)

garch.control(maxiter = 200, trace = TRUE, start = NULL,
  grad = c("analytical","numerical"), abstol = max(1e-20, .Machine$double.eps^2),
  reltol = max(1e-10, .Machine$double.eps^(2/3)), xtol = sqrt(.Machine$double.eps),
  falsetol = 1e2 * .Machine$double.eps, \dots)
}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{order}{a two dimensional integer vector giving the orders of the
    model to fit.  \code{order[2]} corresponds to the ARCH part and
    \code{order[1]} to the GARCH part.}
  \item{series}{name for the series. Defaults to
    \code{deparse(substitute(x))}.}
  \item{control}{a list of control parameters as set up by \code{garch.control}.}
  \item{maxiter}{gives the maximum number of log-likelihood function
    evaluations \code{maxiter} and the maximum number of iterations
    \code{2*maxiter} the optimizer is allowed to compute.}
  \item{trace}{logical. Trace optimizer output?}
  \item{start}{If given this numeric vector is used as the initial estimate
    of the GARCH coefficients.  Default initialization is to set the
    GARCH parameters to slightly positive values and to initialize the
    intercept such that the unconditional variance of the initial GARCH
    is equal to the variance of \code{x}.}
  \item{grad}{character indicating whether analytical gradients or a numerical
    approximation is used for the optimization.}
  \item{abstol}{absolute function convergence tolerance.}
  \item{reltol}{relative function convergence tolerance.}
  \item{xtol}{coefficient-convergence tolerance.}
  \item{falsetol}{false convergence tolerance.}
  \item{\dots}{additional arguments for \code{\link{qr}} when computing
    the asymptotic covariance matrix.}
}
\details{
  \code{garch} uses a Quasi-Newton optimizer to find the maximum
  likelihood estimates of the conditionally normal model.  The first
  max(p, q) values are assumed to be fixed.  The optimizer uses a hessian
  approximation computed from the BFGS update.  Only a Cholesky factor
  of the Hessian approximation is stored.  For more details see Dennis
  et al. (1981), Dennis and Mei (1979), Dennis and More (1977), and
  Goldfarb (1976).  The gradient is either computed analytically or
  using a numerical approximation.
}
\value{
  A list of class \code{"garch"} with the following elements:
  \item{order}{the order of the fitted model.}
  \item{coef}{estimated GARCH coefficients for the fitted model.}
  \item{n.likeli}{the negative log-likelihood function evaluated at the
    coefficient estimates (apart from some constant).}
  \item{n.used}{the number of observations of \code{x}.}
  \item{residuals}{the series of residuals.}
  \item{fitted.values}{the bivariate series of conditional standard
    deviation predictions for \code{x}.}
  \item{series}{the name of the series \code{x}.}
  \item{frequency}{the frequency of the series \code{x}.}
  \item{call}{the call of the \code{garch} function.}
  \item{vcov}{outer product of gradient estimate of the asymptotic-theory
    covariance matrix for the coefficient estimates.}
}
\references{
  A. K. Bera and M. L. Higgins (1993):
  ARCH Models: Properties, Estimation and Testing.
  \emph{J. Economic Surveys} \bold{7} 305--362.
  
  T. Bollerslev (1986):
  Generalized Autoregressive Conditional Heteroscedasticity.
  \emph{Journal of Econometrics} \bold{31}, 307--327.

  R. F. Engle (1982):
  Autoregressive Conditional Heteroscedasticity with Estimates of the
  Variance of United Kingdom Inflation.
  \emph{Econometrica} \bold{50}, 987--1008.

  J. E. Dennis, D. M. Gay, and R. E. Welsch (1981):
  Algorithm 573 --- An Adaptive Nonlinear Least-Squares
  Algorithm.
  \emph{ACM Transactions on Mathematical Software} \bold{7}, 369--383.

  J. E. Dennis and H. H. W. Mei (1979):
  Two New Unconstrained Optimization Algorithms which use Function and
  Gradient Values.
  \emph{J. Optim. Theory Applic.} \bold{28}, 453--482.

  J. E. Dennis and J. J. More (1977):
  Quasi-Newton Methods, Motivation and Theory.
  \emph{SIAM Rev.} \bold{19}, 46--89.

  D. Goldfarb (1976):
  Factorized Variable Metric Methods for Unconstrained Optimization.
  \emph{Math. Comput.} \bold{30}, 796--811.
}
\author{
  A. Trapletti, the whole GARCH part;
  D. M. Gay, the FORTRAN optimizer
}
\seealso{
  \code{\link{summary.garch}} for summarizing GARCH model fits;
  \code{\link{garch-methods}} for further methods.
}
\examples{
n <- 1100
a <- c(0.1, 0.5, 0.2)  # ARCH(2) coefficients
e <- rnorm(n)  
x <- double(n)
x[1:2] <- rnorm(2, sd = sqrt(a[1]/(1.0-a[2]-a[3]))) 
for(i in 3:n)  # Generate ARCH(2) process
{
  x[i] <- e[i]*sqrt(a[1]+a[2]*x[i-1]^2+a[3]*x[i-2]^2)
}
x <- ts(x[101:1100])
x.arch <- garch(x, order = c(0,2))  # Fit ARCH(2) 
summary(x.arch)                     # Diagnostic tests
plot(x.arch)                        

data(EuStockMarkets)  
dax <- diff(log(EuStockMarkets))[,"DAX"]
dax.garch <- garch(dax)  # Fit a GARCH(1,1) to DAX returns
summary(dax.garch)       # ARCH effects are filtered. However, 
plot(dax.garch)          # conditional normality seems to be violated
}
\keyword{ts}
