\name{garch}
\title{Fit GARCH Models to Time Series}
\usage{
garch (x, order = c(1, 1), coef = NULL, itmax = 200, eps = NULL,
       grad = c("analytical","numerical"), series = NULL, trace = TRUE, \dots)
predict (obj, newdata, genuine = FALSE)
coef (obj)
residuals (obj)
fitted (obj)
print (obj, digits = max(3,.Options$digits-3))
summary (obj)
plot (obj)
print.summary (obj, digits = max(3,.Options$digits-3),
               signif.stars = .Options$show.signif.stars, \dots)
}
\alias{garch}
\alias{predict.garch}
\alias{coef.garch}
\alias{residuals.garch}
\alias{fitted.garch}
\alias{print.garch}
\alias{summary.garch}
\alias{plot.garch}
\alias{print.summary.garch}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{order}{a two dimensional integer vector giving the orders of the
    model to fit. \code{order[2]} corresponds to the ARCH part and
    \code{order[1]} to the GARCH part.}
  \item{coef}{If given this numeric vector is used as the initial estimate
    of the GARCH coefficients. Default initialization is to set the
    GARCH parameters to slightly positive values and to initialize the
    intercept such that the unconditional variance of the initial GARCH
    is equal to the variance of \code{x}.}
  \item{itmax}{gives the maximum number of log-likelihood function
    evaluations \code{itmax} and the maximum number of iterations
    \code{2*itmax} the optimizer is allowed to compute.}
  \item{eps}{defines the absolute (\code{max(1e-20,eps^2)}) and
    relative function convergence tolerance
    (\code{max(1e-10,eps^(2/3))}), the coefficient-convergence
    tolerance (\code{sqrt(eps)}), and the false convergence tolerance
    (\code{1e2*eps}). Default value is the machine epsilon, see
    \code{\link{Machine}}.}
  \item{grad}{indicates if the analytical gradient or a numerical
    approximation is used for the optimization.}
  \item{series}{name for the series. Defaults to
    \code{deparse(substitute(x))}.}
  \item{trace}{trace optimizer output?}
  \item{\dots}{additional arguments for \code{\link{qr}} when computing
    the asymptotic standard errors of \code{coef}.}
  \item{obj}{a fit from \code{garch}.}
  \item{newdata}{a numeric vector or time series to compute GARCH
    predictons. Defaults to \code{eval(parse(text=obj$series))}.}
  \item{genuine}{a logical indicating whether a genuine prediction
    should be made, i.e., a prediction for which there is no
    target observation available.}
  \item{digits, signif.stars}{see \code{\link{print.coefmat}}.}
  \item{\dots}{additional arguments for \code{\link{print}}.} 
}
\description{
  Fit a Generalized Autoregressive Conditional Heteroscedastic
  GARCH (p, q) time series model to the data by computing the
  maximum-likelihood estimates of the conditionally normal model.
}
\details{
  \code{garch} uses a Quasi-Newton optimizer to find the
  maximum-likelihood estimates of the conditionally normal model. The
  first max(p,q) values are assumed to be fixed. The optimizer uses a
  hessian approximation computed from the BFGS
  update. Only a Cholesky factor of the Hessian approximation is
  stored. For more details see Dennis et al. (1981), Dennis and Mei
  (1979), Dennis and More (1977), and Goldfarb (1976). The gradient is
  either computed analytically or using a numerical approximation.

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

  \code{print}, \code{summary}, \code{plot}, and \code{print.summary}
  provide the usual generic functions for fitted GARCH
  models. \code{summary} computes the asymptotic standard errors of the
  coefficient estimates from an outer-product approximation of the
  Hessian evaluated at the estimates, see Bollerslev (1986). It furthermore
  tests the residuals for normality and remaining ARCH effects, see
  \code{\link{jarque.bera.test}} and \code{\link{Box.test}}.
  \code{plot} graphically investigates normality and remaining ARCH
  effects for the residuals. 
}
\value{
  For \code{garch} and its methods \code{print} and \code{plot} a list
  of class \code{"garch"} with the following elements:
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
  \item{asy.se.coef}{the asymptotic-theory standard errors of the
    coefficient estimates.}
  
  For \code{predict} a bivariate time series (two-column matrix) of
  predictions. 
  
  For \code{coef}, a numeric vector, for \code{residuals} and
  \code{fitted} a univariate (vector) and a bivariate time series
  (two-column matrix), respectively.

  For \code{summary} and \code{print.summary} a list of
  class \code{"summary.garch"}.
}
\references{
  A. K. Bera and M. L. Higgins (1993): ARCH Models: Properties,
  Estimation and Testing. \emph{J. Economic Surveys} \bold{7} 305-362.
  
  T. Bollerslev (1986): Generalized Autoregressive Conditional 
  Heteroscedasticity. \emph{Journal of Econometrics} \bold{31},
  307-327.

  R. F. Engle (1982): Autoregressive Conditional Heteroscedasticity with
  Estimates of the Variance of United Kingdom
  Inflation. \emph{Econometrica} \bold{50}, 987-1008.

  J. E. Dennis, D. M. Gay, and R. E. Welsch (1981): Algorithm 573 -
  An Adaptive Nonlinear Least-Squares Algorithm. \emph{ACM Trans. 
  Math. Software} \bold{7}, 369-383.

  J. E. Dennis and H. H. W. Mei (1979): Two New Unconstrained Optimization
  Algorithms which use Function and Gradient
  Values. \emph{J. Optim. Theory Applic.} \bold{28}, 453-482.

  J. E. Dennis and J. J. More (1977): Quasi-Newton Methods, Motivation
  and Theory. \emph{SIAM Rev.} \bold{19}, 46-89.

  D. Goldfarb (1976): Factorized Variable Metric Methods for
  Unconstrained Optimization. \emph{Math. Comput.} \bold{30}, 796-811.
}
\author{
  A. Trapletti, the whole GARCH part, D. M. Gay, the fortran optimizer
}
\examples{
n <- 1100
a <- c (0.1, 0.5, 0.2)  # ARCH(2) coefficients
e <- rnorm (n)  
x <- double (n)
x[1:2] <- rnorm (2, sd = a[1]/(1.0-a[2]-a[3])) 
for (i in 3:n)  # Generate ARCH(2) process
{
  x[i] <- e[i]*sqrt(a[1]+a[2]*x[i-1]^2+a[3]*x[i-2]^2)
}
x <- ts(x[101:1100])
x.arch <- garch (x, order = c(0,2))  # Fit ARCH(2) 
summary (x.arch)                     # Diagnostic tests
plot (x.arch)                        

data (EuStockMarkets)  
dax <- diff(log(EuStockMarkets))[,"DAX"]
dax.garch <- garch (dax)  # Fit a GARCH(1,1) to DAX returns
summary (dax.garch)       # ARCH effects are filtered. However, 
plot (dax.garch)          # conditional normality seems to be violated
}
\keyword{ts}


