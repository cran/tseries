\name{bootstrap}
\alias{bootstrap}
\alias{print.resample.statistic}
\title{Bootstrap for General Stationary Data}
\description{
  \code{bootstrap} generates bootstrap samples for general stationary
  data and computes the bootstrap estimate of standard error and bias
  if a statistic is given.
}
\usage{
bootstrap(x, nb = 1, statistic = NULL, b = NULL,
          type = c("stationary","block"), \dots)
}
\arguments{
  \item{x}{a numeric vector or time series giving the original data.}
  \item{nb}{the number of bootstrap series to compute.}
  \item{statistic}{a function which when applied to a time series
    returns a vector containing the statistic(s) of interest.}
  \item{b}{if \code{type} is \code{"stationary"}, then \code{b} is the
    mean block length. If \code{type} is \code{"block"}, then \code{b}
    is the fixed block length.}
  \item{type}{the type of bootstrap to generate the simulated time
    series. The possible input values are \code{"stationary"}
    (stationary bootstrap with mean block length \code{b}) and
    \code{"block"} (blockwise bootstrap with block length
    \code{b}). Default to \code{"stationary"}.} 
  \item{\dots}{additional arguments for \code{statistic} which are
    passed unchanged each time \code{statistic} is called.}
}
\details{
  If \code{type} is \code{"stationary"}, then the stationary
  bootstrap scheme with mean block length \code{b} according to Politis
  and Romano (1994) is used to generate the simulated series. For
  \code{type} equals \code{"block"}, the blockwise bootstrap with block
  length \code{b} according to Kuensch (1989) is used.

  Note, that this bootstrap procedure does not implement the block of
  blocks bootstrap. Hence, it may only be used for symmetric statistics,
  i.e., for statistics which are invariant under permutations of
  \code{x}.
  
  For consistency, the (mean) block length \code{b} should grow with
  \code{n} at an appropriate rate. If \code{b} is not given, then a
  default growth rate of \code{const * n^(1/3)} is used (\code{n} is the
  number of observations in \code{x}). This rate is "optimal" under
  certain conditions (see the references for more details). However, in
  general the growth rate depends on the specific properties of the
  process which has generated \code{x}. Therefore, a default value for
  \code{const} has been determined by a Monte Carlo simulation using a
  Gaussian AR(1) process (AR(1)-parameter of 0.5, 500
  observations). \code{const} has been chosen such that the mean square
  error for the bootstrap estimate of the variance of the empirical mean
  is minimized. 

  Missing values are not allowed.

  There is a special print method for objects of class
  \code{"resample.statistic"} which by default uses
  \code{max(3, getOption("digits") - 3)} digits to format real numbers.
}
\value{
  If \code{statistic} is \code{NULL}, then it returns a matrix or time
  series with \code{nb} columns and \code{length(x)} rows containing the
  bootstrap data. Each column contains one bootstrap sample.

  If \code{statistic} is given, then a list of class
  \code{"resample.statistic"} with the following elements is returned:
  \item{statistic}{the results of applying \code{statistic} to each of
    the simulated time series.}
  \item{orig.statistic}{the results of applying \code{statistic} to the
    original series.}
  \item{bias}{the bootstrap estimate of the bias of \code{statistic}.}
  \item{se}{the bootstrap estimate of the standard error of \code{statistic}.}
  \item{call}{the original call of \code{bootstrap}.}
}
\references{
  H. R. Kuensch (1989):
  The Jackknife and the Bootstrap for General Stationary Observations.
  \emph{The Annals of Statistics} \bold{17}, 1217--1241.

  D. N. Politis and J. P. Romano (1994):
  The Stationary Bootstrap.
  \emph{Journal of the American Statistical Association} \bold{89},
  1303--1313.
}
\author{A. Trapletti}
\seealso{
  \code{\link{sample}},
  \code{\link{surrogate}}
}
\examples{
n <- 500  # Generate AR(1) process
a <- 0.6
e <- rnorm(n+100)  
x <- double(n+100)
x[1] <- rnorm(1)
for(i in 2:(n+100)) {
  x[i] <- a * x[i-1] + e[i]
}
x <- ts(x[-(1:100)])

bootstrap(x, nb=500, statistic=mean)  # ok 
bootstrap(x, nb=500, statistic=mean, type="block", b=1)  # wrong here

# std. error for mean should approximately be
sqrt(1/(n*(1-a)^2))
}
\keyword{ts}
