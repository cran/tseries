# Copyright (C) 1997-1999  Adrian Trapletti
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the Free
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#
# Descriptive time series analysis in the time domain
#


acf <- function (x, lag = length(x)-1, correlation = TRUE, pl = TRUE, ...)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  if (lag >= length(x)) stop ("Number of lags should be less than number of data points")
  if (lag < 1) stop ("Number of lags is not positive")
  x <- scale (x, scale = FALSE)
  n <- length (x)
  x <- c (x, rep (0, nextn(n+lag)-n))
  a <- fft (x)
  ac <- Re (fft (a*Conj(a),inv=TRUE)/length(x))
  acov <- ac[1:(lag+1)]/n
  lg <- 0:lag
  if (correlation)
  {
    wb <- rep (1.96/sqrt(n), lag+1)
    acf <- acov/acov[1]
    code <- c ("", "*")
    ysgnf <- code[as.integer(((wb<acf)|(-wb>acf)) & c(F,rep(T,length(lg)-1)))+1]
    acf <- tsparam (x=lg, y=acf, xlab="lag", ylab="autocorrelation function", ybnd=list(wb,-wb),
                    ysgnf=ysgnf)
  }
  else 
    acf <- tsparam (x=lg, y=acov, xlab="lag", ylab="autocovariance function")
  if (pl) plot (acf, ...)
  return (acf)
}

ccf <- function (x, y, lag = length(x)-1, correlation = TRUE, pl = TRUE, ...)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  if (!is.vector(y) & !is.univariate.ts(y))
    stop ("x is not a vector or univariate time series")
  if (length(x) != length(y)) stop ("x and y have not the same length")
  if ((is.ts(x) & !is.ts(y)) | (!is.ts(x) & is.ts(y)))
    stop ("x and y must have the same type")
  if (is.ts(x))
    if (!all(tsp(x) == tsp(y)))
      stop ("x and y must be aligned for crosscorrelation analysis")
  if (lag >= length(x)) stop ("Number of lags should be less than number of data points")
  if (lag < 1) stop ("Number of lags is not positive")
  x <- scale (x, scale = FALSE)
  y <- scale (y, scale = FALSE)
  n <- length (x)
  x <- c (x, rep (0, nextn(n+lag)-n))
  y <- c (y, rep (0, nextn(n+lag)-n))
  a <- fft (x)
  b <- fft (y)
  cc <- Re (fft (a*Conj(b),inv=TRUE)/length(x))
  ccov <- c (cc[(length(cc)-lag+1):length(cc)]/n, cc[1:(lag+1)]/n)
  lg <- -lag:lag
  if (correlation)
  {
    ccx <- Re (fft (a*Conj(a),inv=TRUE)/length(x))
    ccy <- Re (fft (b*Conj(b),inv=TRUE)/length(x))
    ccovx1 <- ccx[1]/n
    ccovy1 <- ccy[1]/n
    wb <- rep (1.96/sqrt(n), 2*lag+1)
    ccf <- ccov/sqrt(ccovx1*ccovy1)
    code <- c ("", "*")
    ysgnf <- code[as.integer(((wb<ccf)|(-wb>ccf)) & c(F,rep(T,length(lg)-1)))+1]
    ccf <- tsparam (x=lg, y=ccf, xlab="lag", ylab="crosscorrelation function", ybnd=list(wb,-wb),
                    ysgnf=ysgnf)
  }
  else 
    ccf <- tsparam (x=lg, y=ccov, xlab="lag", ylab="crosscovariance function")
  if (pl) plot (ccf, ...)
  return (ccf)
}

amif <- function (x, lag = length(x)-1, maxbit = 20, confidence = 0.2, nsurr = 20,
                  fft = FALSE, amplitude = FALSE, normalized = TRUE, trace = FALSE, pl = TRUE, ...)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  if (lag >= length(x)) stop ("Number of lags should be less than number of data points")
  if (lag < 1) stop ("Number of lags is not positive")
  if ((maxbit < 1) | (maxbit > 25)) stop ("maxbit out of range")
  if ((confidence < 0.01) | (confidence > 0.99)) stop ("confidence out of range")
  if (nsurr < 1) stop ("nsurr is not positive")
  lx <- length (x)
  inf <- double (lag+1)
  res <- .C ("R_amif", as.vector(x,mode="double"), as.integer(lx), inf=as.vector(inf),
             as.integer(lag), as.integer(maxbit), as.double(confidence),
             as.integer(normalized), as.integer(trace))
  cor <- res$inf
  surrsam <- surrogate (x, ns=nsurr, fft=fft, amplitude=amplitude)
  surrwb <- matrix (0, nrow = nsurr, ncol = lag+1)
  for (i in 1:nsurr)
  {
    res <- .C ("R_amif", as.vector(surrsam[,i],mode="double"), as.integer(lx), inf=as.vector(inf),
               as.integer(lag), as.integer(maxbit), as.double(confidence),
               as.integer(normalized), as.integer(trace))
    surrwb[i,] <- res$inf
  }
  wb <- apply (surrwb, 2, quantile, 0.95) 
  lg <- 0:lag
  code <- c ("", "*")
  ysgnf <- code[as.integer((wb<cor) & c(F,rep(T,length(lg)-1)))+1]
  if (normalized) ylab <- "normalized auto mutual information function"
  else ylab <- "auto mutual information function"
  amif <- tsparam (x=lg, y=cor, xlab="lag", ylab=ylab, ybnd=list(wb),
                   ysgnf=ysgnf)
  if (pl) plot (amif, ...)
  return (amif)
}

pacf <- function (x, lag = length(x) - 1, pl = TRUE, ...)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  if (lag >= length(x)) stop ("Number of lags should be less than number of data points")
  if (lag < 1) stop ("Number of lags is not positive")
  n <- length (x)
  cr <- acf (x, lag=lag, correlation=TRUE, pl=FALSE)
  pacf <- matrix (0, lag, lag)
  cor <- cr$y[2:(lag+1)]
  res <- .C ("R_Durbin", as.vector(cor,mode="double"),
             pacf=as.matrix(pacf), as.integer(lag))
  pacf <- diag(res$pacf)
  pacf <- c(1,pacf)
  lg <- 0:lag
  wb <- rep (1.96/sqrt(n), lag+1)
  code <- c ("", "*")
  ysgnf <- code[as.integer(((wb<pacf)|(-wb>pacf)) & c(F,rep(T,length(lg)-1)))+1]
  pacf <- tsparam (x=lg, y=pacf, xlab="lag", ylab="partial autocorrelation function",
                   ybnd=list(wb,-wb), ysgnf=ysgnf)
  if (pl) plot (pacf, ...)
  return (pacf)
}
