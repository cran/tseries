\name{tsbootstrap}
\alias{tsbootstrap}
\alias{print.resample.statistic}
\title{Bootstrap for General Stationary Data}
\description{
  \code{tsbootstrap} generates bootstrap samples for general stationary
  data and computes the bootstrap estimate of standard error and bias
  if a statistic is given.
}
\usage{
tsbootstrap(x, nb = 1, statistic = NULL, m = 1, b = NULL,
            type = c("stationary","block"), \dots)
}
\arguments{
  \item{x}{a numeric vector or time series giving the original data.}
  \item{nb}{the number of bootstrap series to compute.}
  \item{statistic}{a function which when applied to a time series
    returns a vector containing the statistic(s) of interest.}
  \item{m}{the length of the basic blocks in the block of blocks
    bootstrap.} 
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
  and Romano (1994) is computed. For \code{type} equals \code{"block"},
  the blockwise bootstrap with block length \code{b} according to
  Kuensch (1989) is used. 

  If \code{m > 1}, then the block of blocks bootstrap is computed
  (see Kuensch, 1989). The basic sampling scheme is the same as for 
  the case \code{m = 1}, except that the bootstrap is applied to a series
  \code{y} containing blocks of length \code{m}, where each block of \code{y} is
  defined as \eqn{y[t] = (x[t], \dots, x[t-m+1])}. Therefore, for the block
  of blocks bootstrap the first argument of \code{statistic} is given by
  a \code{n x m} matrix \code{yb}, where each row of \code{yb} contains one
  bootstrapped basic block observation \eqn{y[t]} (\code{n} is the number of
  observations in \code{x}). 

  Note, that for statistics which are functions of the empirical
  \code{m}-dimensional marginal (\code{m > 1}) only this procedure
  yields asymptotically valid bootstrap estimates. The 
  case \code{m = 1} may only be used for symmetric statistics (i.e., for
  statistics which are invariant under permutations of \code{x}).
  \code{\link[boot]{tsboot}} does not implement the block of blocks
  bootstrap, and, therefore, the first example in \code{\link[boot]{tsboot}}
  yields inconsistent estimates.
  
  For consistency, the (mean) block length \code{b} should grow with
  \code{n} at an appropriate rate. If \code{b} is not given, then a
  default growth rate of \code{const * n^(1/3)} is used. This rate is
  "optimal" under certain conditions (see the references for more
  details). However, in general the growth rate depends on the specific
  properties of the data generation process. A default value for
  \code{const} has been determined by a Monte Carlo simulation using a
  Gaussian AR(1) process (AR(1)-parameter of 0.5, 500
  observations). \code{const} has been chosen such that the mean square
  error for the bootstrap estimate of the variance of the empirical mean
  is minimized.  

  Note, that the computationally intensive parts are fully implemented
  in \code{C} which makes \code{tsbootstrap} about 10 to 30 times faster
  than \code{\link[boot]{tsboot}}.  
  
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
  \item{call}{the original call of \code{tsbootstrap}.}
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
  \code{\link{surrogate}},
  \code{\link[boot]{tsboot}}
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

tsbootstrap(x, nb=500, statistic=mean)

# Asymptotic formula for the std. error of the mean
sqrt(1/(n*(1-a)^2))

acflag1 <- function(x)
{
  xo <- c(x[,1], x[1,2])
  xm <- mean(xo)
  return(mean((x[,1]-xm)*(x[,2]-xm))/mean((xo-xm)^2))
}

tsbootstrap(x, nb=500, statistic=acflag1, m=2)

# Asymptotic formula for the std. error of the acf at lag one
sqrt(((1+a^2)-2*a^2)/n)
}
\keyword{ts}
