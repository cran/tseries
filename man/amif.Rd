\name{amif}
\title{Auto Mutual Information Function}
\usage{
amif (x, lag.max = NULL, maxbit = 20, confidence = 0.2, ci = 0.95, nsurr = 20, 
      fft = FALSE, amplitude = FALSE, normalized = TRUE, trace = FALSE,
      plot = TRUE, na.action = na.fail, \dots)
plot (obj, ci.col = "blue", ...)
}
\alias{amif}
\alias{plot.amif}
\arguments{
  \item{x}{a numeric vector or time series.}
  \item{lag.max}{a scalar lag parameter.}
  \item{maxbit}{the maximum resolution in bit to which the precision of
    the data will be limited; also limits the maximum partition depth. It 
    cannot exceed an overall maximum of 25 bit.}
  \item{confidence}{the confidence level for the chi-square-test 
    which tests the point distribution in substructures of the
    probability area upon uniformity. A substructure is assumed if the
    deviation from a uniform point distribution occurs with a probability
    less than the significance niveau. Possible settings are:
    0.01, 0.02, 0.05, 0.1, 0.2, \dots, 0.9, 0.95, 0.98, 0.99.}
  \item{ci}{coverage probability for confidence bound. Plotting of
    the confidence bound is suppressed if \code{ci} is zero or negative.}
  \item{nsurr}{the number of surrogate samples to compute the confidence
    bound. To obtain "good" confidence bounds, set \code{nsurr} to at least
    a value of \code{500}.}
  \item{fft}{a logical indicating whether phase randomized surrogate
    data is generated.}
  \item{amplitude}{a logical indicating whether amplitude-adjusted
    surrogate data is computed.}
  \item{normalized}{a logical indicating whether the normalized auto 
    mutual information function is computed.}
  \item{trace}{a logical indicating whether additional output from the
    computation is traced.}
  \item{plot}{a logical indicating whether the AMIF is plotted.}
  \item{na.action}{function to be called to handle missing values.}
  \item{\dots}{additional arguments to \code{\link{plot.amif}}.}
  \item{obj}{an object of class \code{"amif"}.}
  \item{ci.col}{color to plot the confidence bound.}
  \item{\dots}{additional arguments to \code{\link{plot.acf}}.}
}
\description{
  Computes and plots the sample (normalized) auto mutual information function
  of \code{x} up to lag \code{lag.max}. If \code{plot} is \code{TRUE}, then
  the AMIF is plotted. Also the confidence bound for the null of
  surrogate data is computed by Monte-Carlo simulation percentiles
  and plotted. If any temporal dependence is
  under question, then the null is an i.i.d. series and scrambled
  surrogates can be used (the default). In the case of testing
  non-linearity against linearity, setting the switch \code{fft} is
  appropriate. It generates surrogates with the same spectrum as
  \code{x}. If the switch \code{amplitude} is set in addition, then
  surrogates \code{xs} with the following properties are used: First, \code{xs}
  has the same histogram as \code{x}. Second, \code{G(xs)} has the 
  same Fourier spectra as \code{G(x)}, where \code{G(.)} is the
  transformation from the histogram of \code{x} to a Gaussian
  distribution.

  The simulations are computed with the actual values of 
  \code{maxbit} and \code{confidence}. 
}
\details{
  To compute the AMIF, an implementation of the algorithm of Fraser and
  Swinney (1986) is used.
}
\value{
  An object of class \code{"amif","acf"}.
}
\references{
  C. Granger and J. L. Lin (1994): Using the Mutual Information
  Coefficient to Identify Lags in Non-Linear Models. \emph{Journal of
  Time Series Analysis} \bold{15}, pp. 371-384.

  A. M. Fraser and H. L. Swinney (1986): Independent Coordinates for
  Strange Attractors from Mutual Information. \emph{Physical Review A}
  \bold{33}, pp. 1134-1140.
}
\author{Fraser and Swinney (1986) algorithm: 
  The group of F. W. Schneider at the University of Wuerzburg, Germany:
  Authors T. M. Kruel, Institut fuer Physikalische Chemie,
  University of Wuerzburg, Germany, and K. Krischer, 
  Fritz-Haber-Institut der Max-Planck-Gesellschaft, Berlin, Germany.

  Port to R, implementation of AMIF, and scramble device: A. Trapletti.
}
\seealso{
  \code{\link{acf}},
  \code{\link{surrogate}} 
}
\examples{
n <- 1000  # Length of simulated time series

e <- rnorm (n)  # Generate ARCH(1) process
x <- double (n)
x[1] <- rnorm (1)
for (i in 2:n)
{
  x[i] <- e[i]*sqrt(0.1+0.4*x[i-1]^2)
}
x <- ts(x)
plot (x)

# Each test takes about 3 sec on a Pentium II 300MHz

amif (x, lag.max=5)  # i.i.d. vs. any dependence 
amif (x, lag.max=5, fft=T)  # linear vs. non-linear 
amif (x, lag.max=5, fft=T, amplitude=T)

e <- rnorm (n)  # Generate AR(1) process
x <- double (n)
x[1] <- rnorm (1)
for (i in 2:n)
{
  x[i] <- 0.4*x[i-1]+e[i]
}
x <- ts(x)
plot (x)

amif (x, lag.max=5)  # i.i.d. vs. any dependence 
amif (x, lag.max=5, fft=T)  # linear vs. non-linear 
amif (x, lag.max=5, fft=T, amplitude=T)
}
\keyword{ts}

