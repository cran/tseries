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
# Frequency domain time series analysis
#


spectrum <- function (x, k = fejer.kernel(nextn(length(x))%/%2-1,length(x)/10), pl = TRUE, ...)
{
  if (!is.vector(x) & !is.univariate.ts(x)) 
    stop("x is not a vector or univariate time series")
  if (!is.kernel(k))
    stop ("k is not a kernel")
  x <- scale (x, scale = F)  # Compute periodogram
  n <- length (x)
  npad <- nextn (n)
  x <- c (x, rep(0,npad-n))  
  per <- (Mod(fft(x))^2)/(2*pi*n)  # Divide by n to preserve the area under the periodogram
  per[1] <- 0.0
  per <- apply.kernel (per, k, circular=T)  # smooth periodogram with kernel k    
  freq <- seq (0, 0.5, length=(npad%/%2+1))
  df <- df.kernel(k)
  spectrum <- tsparam(x=freq, y=per[1:(npad%/%2+1)], xlab="frequency", ylab="spectrum",
                      lab=paste(attr(k,"name"),
                        ", bw ", formatC(df/(2*n),digits=3,format="f"),
                        ", df ", formatC(df,format="d")))
  if (pl) plot (spectrum, ...)
  return (spectrum)
}

cumulative.periodogram <- function (x, pl = TRUE, ...)
{
  if (!is.vector(x) & !is.univariate.ts(x)) 
    stop("x is not a vector or univariate time series")
  s <- spectrum (x, daniell.kernel(0), pl=F)
  cum <- cumsum(s$y)/sum(s$y)
  q <- floor((length(x)-1)/2)
  k1 <- (2*s$x*q-1)/(q-1)
  k2 <- 1.36/sqrt(q-1)
  wb1 <- k1+k2
  wb2 <- k1-k2
  code <- c ("", "*")
  ysgnf <- code[as.integer((wb1<cum)|(wb2>cum))+1]
  cmpgram <- tsparam(x=s$x, y=cum, xlab="frequency",
                     ylab="cumulative periodogram",ybnd=list(wb1,wb2), ysgnf=ysgnf)
  if (pl) plot (cmpgram, ...)
  return (cmpgram)
}

cross.spectrum <- function (x, y, k = fejer.kernel(nextn(length(x))%/%2-1,length(x)/10))
{
  if (!is.vector(x) & !is.univariate.ts(x)) 
    stop ("x is not a vector or univariate time series")
  if (!is.vector(y) & !is.univariate.ts(y)) 
    stop ("y is not a vector or univariate time series")
  if (length(x) != length(y)) 
    stop("x and y have not the same length")
  if ((is.ts(x) & !is.ts(y)) | (!is.ts(x) & is.ts(y)))
    stop ("x and y must have the same type")
  if (is.ts(x))
    if (!all(tsp(x) == tsp(y)))
      stop ("x and y must be aligned for cross-spectral analysis")
  if (!is.kernel(k))
    stop ("k is not a kernel")
  x <- scale (x, scale = F)  # Compute cross periodogram
  y <- scale (y, scale = F)
  n <- length (x)
  npad <- nextn (n)
  x <- c (x, rep(0,npad-n))
  y <- c (y, rep(0,npad-n))
  per <- fft(x)*Conj(fft(y))/(2*pi*n)  
  per[1] <- 0.0
  per <- apply.kernel (per, k, circular=T)  # smooth periodogram with kernel k    
  freq <- seq (0, 0.5, length=(npad%/%2+1))
  df <- df.kernel(k)
  spectrum <- tsparam(x=freq, y=per[1:(npad%/%2+1)], xlab="frequency", ylab="cross spectrum",
                      lab=paste(attr(k,"name"),
                        ", bw ", formatC(df/(2*n),digits=3,format="f"),
                        ", df ", formatC(df,format="d")))
  return (spectrum)
}

co.spectrum <- function (x, y, k = fejer.kernel(nextn(length(x))%/%2-1,length(x)/10),
                         pl = TRUE, ...)
{
  s <- cross.spectrum(x,y,k)
  s$y <- Re(s$y)
  attr(s, "ylab") <- "cospectrum"
  if (pl) plot (s, ...)
  return (s)
}

quadrature.spectrum <- function (x, y, k = fejer.kernel(nextn(length(x))%/%2-1,length(x)/10),
                                 pl = TRUE, ...)
{
  s <- cross.spectrum(x,y,k)
  s$y <- -Im(s$y)
  attr(s, "ylab") <- "quadrature spectrum"
  if (pl) plot (s, ...)
  return (s)
}

absolute.coherency <- function (x, y, k = fejer.kernel(nextn(length(x))%/%2-1,length(x)/10),
                                pl = TRUE, ...)
{
  sx <- spectrum(x,k,pl=F)
  sy <- spectrum(y,k,pl=F)
  sxy <- cross.spectrum(x,y,k)
  spec <- Mod(sxy$y)/sqrt(sx$y*sy$y)
  df <- df.kernel(k)
  if (k$m < 1)
    ybnd <- NULL
  else
  {
    b <- 1.96*(1-spec^2)/(sqrt(df))
    wb1 <- spec+b
    wb2 <- spec-b
    ybnd <- list(wb1,wb2)
  }
  abscoh <- tsparam (x=sxy$x, y=spec, xlab="frequency", ylab="absolute coherency",
                     ybnd=ybnd,
                     lab=paste(attr(k,"name"),
                       ", bw ", formatC(df/(2*length(x)),digits=3,format="f"),
                       ", df ", formatC(df,format="d")))
  if (pl) plot (abscoh, ...)
  return (abscoh)
}

phase.spectrum <- function (x, y, k = fejer.kernel(nextn(length(x))%/%2-1,length(x)/10),
                            pl = TRUE, ...)
{
  sx <- spectrum(x,k,pl=F)
  sy <- spectrum(y,k,pl=F)
  sxy <- cross.spectrum(x,y,k)
  phspec <- Arg(sxy$y)
  abscoh <- Mod(sxy$y)/sqrt(sx$y*sy$y)
  df <- df.kernel(k)
  if (k$m < 1)
    ybnd <- NULL
  else
  {
    b <- 1.96*sqrt(1/abscoh^2-1)/(sqrt(df))
    wb1 <- phspec+b
    wb2 <- phspec-b
    ybnd <- list(wb1,wb2)
  }
  phspec <- tsparam (x=sxy$x, y=phspec, xlab="frequency", ylab="phase spectrum",
                     ybnd=ybnd,
                     lab=paste(attr(k,"name"),
                       ", bw ", formatC(df/(2*length(x)),digits=3,format="f"),
                       ", df ", formatC(df,format="d")))
  if (pl) plot (phspec, ...)
  return (phspec)
}

amplitude.spectrum <- function (x, y, k = fejer.kernel(nextn(length(x))%/%2-1,length(x)/10),
                                pl = TRUE, ...)
{
  sx <- spectrum(x,k,pl=F)
  sy <- spectrum(y,k,pl=F)
  sxy <- cross.spectrum(x,y,k)
  ampspec <- Mod(sxy$y)
  abscoh <- ampspec/sqrt(sx$y*sy$y)
  df <- df.kernel(k)
  if (k$m < 1)
    ybnd <- NULL
  else
  {
    b <- 1.96*ampspec*sqrt(1/abscoh^2+1)/(sqrt(df))
    wb1 <- ampspec+b
    wb2 <- ampspec-b
    ybnd <- list(wb1,wb2)
  }  
  ampspec <- tsparam (x=sxy$x, y=ampspec, xlab="frequency", ylab="amplitude spectrum",
                      ybnd=ybnd,
                      lab=paste(attr(k,"name"),
                        ", bw ", formatC(df/(2*length(x)),digits=3,format="f"),
                        ", df ", formatC(df,format="d")))
  if (pl) plot (ampspec, ...)
  return (ampspec)
}




