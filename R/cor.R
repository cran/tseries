## Copyright (C) 1997-1999  Adrian Trapletti
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

##
## Descriptive time series analysis in the time domain
##


amif <- function (x, lag.max = NULL, maxbit = 20, confidence = 0.2, ci = 0.95, nsurr = 20,
                  fft = FALSE, amplitude = FALSE, normalized = TRUE, trace = FALSE,
                  plot = TRUE, ...)
{
  if (NCOL(x) > 1) stop ("x is not a vector or univariate time series")
  if ((maxbit < 1) | (maxbit > 25)) stop ("maxbit out of range")
  if ((confidence < 0.01) | (confidence > 0.99)) stop ("confidence out of range")
  if (nsurr < 1) stop ("nsurr is not positive")
  if (ci >= 1) stop ("ci out of range")
  series <- deparse (substitute(x))
  x <- as.ts(x)
  x.freq <- frequency(x)
  x <- as.matrix(x)
  if (any(is.na(x))) stop ("NAs in x")
  sampleT <- nrow(x)
  if (is.null(lag.max))
    lag.max <- floor (10*(log10(sampleT)-log10(1)))
  lag.max <- min (lag.max, sampleT-1)
  if (lag.max < 1) 
    stop ("lag.max must be at least 1")
  lag <- matrix (1, 1, 1)
  lag[lower.tri(lag)] <- -1
  inf <- double (lag.max+1)
  cor <- array(.C ("R_amif", as.vector(x,mode="double"), as.integer(sampleT), inf=as.vector(inf),
                   as.integer(lag.max), as.integer(maxbit), as.double(confidence),
                   as.integer(normalized), as.integer(trace), PACKAGE="tseries")$inf,
               c(lag.max + 1, 1, 1))
  if (ci > 0)
  {
    surrsam <- surrogate (x, ns=nsurr, fft=fft, amplitude=amplitude)
    surrwb <- matrix (0, nrow = nsurr, ncol = lag.max+1)
    for (i in 1:nsurr)
      surrwb[i,] <- .C ("R_amif", as.vector(surrsam[,i],mode="double"), as.integer(sampleT),
                        inf=as.vector(inf), as.integer(lag.max), as.integer(maxbit),
                        as.double(confidence), as.integer(normalized),
                        as.integer(trace), PACKAGE="tseries")$inf
    wb <- apply (surrwb, 2, quantile, ci)
  }
  else
    wb <- NULL
  lag <- outer(0:lag.max, lag/x.freq)
  amif <- structure(.Data = list(acf = cor, type = "covariance", n.used = sampleT, lag = lag,
                      series = series, snames = colnames(x), clim = wb, normalized = normalized),
                    class = c("amif","acf"))
  if (plot)
  {
    plot.amif (amif, ...)
    return (invisible(amif))
  }
  else return (amif)  
}

plot.amif <- function (object, ci.col = "blue", ylab = "AMIF", ...)
{
  if (!inherits(object, "amif")) stop ("method is only for amif objects")
  plot.acf (x = object, ylab = ylab, ...)
  if (!is.null(object$clim))
    lines(object$lag[, 1, 1], object$clim, col = ci.col, lty = 2)
}
