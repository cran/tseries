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
# Various time series related routines
#


intgrt <- function (obj, ...) { UseMethod("intgrt") }

embed <- function (x, dimension = 1) 
{
  if (is.vector(x) | is.univariate.ts(x))
  {
    n <- length (x)
    if ((dimension < 1) | (dimension > n)) stop ("wrong embedding dimension")
    y <- matrix (0.0, n-dimension+1, dimension)
    res <- .C("R_embed_vec", as.vector(x,mode="double"), y=as.matrix(y),
              as.integer(dimension),as.integer(n))
    y <- res$y
  }
  else if (is.matrix(x))
  {
    n <- dim(x)[1]
    m <- dim(x)[2]
    if ((dimension < 1) | (dimension > n)) stop ("wrong embedding dimension")
    y <- matrix (0.0, n-dimension+1, dimension*m)
    for (i in (1:m))
      y[,seq(i,by=m,length=dimension)] <- embed (as.vector(x[,i]), dimension)
  }
  else
    stop ("x is not a vector or matrix")
  return (y)
}

intgrt.vec <- function (x, lag = 1, differences = 1, xi = rep(0.0,lag*differences))
{
  if (!is.vector(x)) stop ("x is not a vector")
  if (lag < 1 | differences < 1) stop ("Bad value for lag or differences")
  if (length(xi) != lag*differences) stop ("xi has not the right length")
  if (differences == 1)
  {
    n <- length(x)
    x <- as.vector(x,mode="double")
    y <- as.vector(numeric(n+lag))
    xi <- as.vector(xi,mode="double")
    for (i in 1:lag)
      y[i] <- xi[i]
    res <- .C ("R_intgrt_vec", x, y=y, as.integer(lag), as.integer(n))
    y <- res$y
  }
  else
    y <- intgrt.vec (intgrt.vec (x, lag, differences-1, xi[(lag+1):(lag*differences)]),
                     lag, 1, xi[1:lag])
  return (y)
}

intgrt.default <- function (x, lag = 1, differences = 1,
                            xi = rep(0.0,lag*differences*dim(as.matrix(x))[2]))
{
  if (is.matrix(x))
  {
    n <- dim(x)[1]
    m <- dim(x)[2]
    y <- matrix (0, nr=n+lag*differences, nc=m)
    dim(xi) <- c(lag*differences,m)
    for (i in 1:m)
      y[,i] <- intgrt.vec (as.vector(x[,i]), lag, differences, as.vector(xi[,i]))
  }
  else if (is.vector(x))
    y <- intgrt.vec (x, lag, differences, xi)
  else
    stop ("x is not a vector or matrix")
  return (y)
}

intgrt.ts <- function (x, lag = 1, differences = 1,
                       xi = rep(0.0,lag*differences*dim(as.matrix(x))[2]))
{
  if (is.univariate.ts(x))
    y <- intgrt.default (as.vector(x), lag, differences, xi)
  else
    y <- intgrt.default (as.matrix(x), lag, differences, xi)
  y <- ts (y, frequency = frequency(x), start = start(x))
  return (y)
}

is.univariate.ts <- function (obj)
{
  return (is.ts(obj) & is.null(dim(obj)))
}

read.ts <- function (file, header = FALSE, sep = "", skip = 0, ...)
{
  x <- read.matrix (file, header = header, sep = sep, skip = skip)
  x <- ts (x, ...)
  return (x)
}

fftsurr <- function (x)
  # This is algorithm 1, p. 183 from "Theiler et al. (1992): Using
  # Surrogate Data to Detect Nonlinearity in Time Series, in Nonlinear
  # Modelling and Forecasting, Editors Casdagli & Eubank, Santa Fe Institute,
  # Addison Wesley". Note that Step 7. and 8. are only for t = 2,...,N.
{
  z <- fft (x);
  zz <- z*exp(1i*runif(z, max=2*pi));  
  re <- Re(zz[2:length(zz)]+zz[length(zz):2])/2
  im <- Im(zz[2:length(zz)]-zz[length(zz):2])/2
  zzz1 <- Re(zz[1]+zz[1])/2+1i*Im(zz[1]-zz[1])/2 
  zzz <- c(zzz1,re+1i*im)
  return (Re(fft(zzz, inverse=T)))
}

ampsurr <- function (x)
  # This is algorithm 2, pp. 183, 184 from "Theiler et al. (1992): Using
  # Surrogate Data to Detect Nonlinearity in Time Series, in Nonlinear
  # Modelling and Forecasting, Editors Casdagli & Eubank, Santa Fe Institute,
  # Addison Wesley". 
{
  sx <- sort(x)
  rx <- rank(x)
  g <- rnorm(x)
  sg <- sort(g)
  y <- sg[rx]
  yy <- fftsurr(y)
  ryy <- rank(yy)
  return (sx[ryy])
}

surrogate <- function (x, ns = 1, fft = FALSE, amplitude = FALSE)
{
  if (!is.vector(x) & !is.univariate.ts(x))
    stop ("x is not a vector or univariate time series")
  if (ns < 1) stop ("ns is not positive")
  n <- length(x)
  surrogate <- matrix (x, nrow=n, ncol=ns)
  if (fft)
  {
    if (amplitude)
      surrogate <- apply(surrogate, 2, ampsurr)
    else
      surrogate <- apply(surrogate, 2, fftsurr)
  }
  else
    surrogate <- apply(surrogate, 2, sample, replace=F)
  return (surrogate)
}

quadmap <- function (xi = 0.0, a = 4.0, n = 1000)
{
  if (n < 1) stop ("n is not positive")
  if ((xi < 0) | (xi > 1)) stop ("xi is not in [0,1]")
  if ((a < 0) | (xi > 4)) stop ("a is not in [0,4]")
  x <- double(n)
  res <- .C ("R_quad_map", x=as.vector(x), as.double(xi), as.double(a), as.integer(n))
  return (ts(res$x))
}

read.matrix <- function (file, header = FALSE, sep = "", skip = 0)
{
  row.lens <- count.fields (file, sep = sep, skip = skip)
  if (any (row.lens != row.lens[1])) 
    stop ("number of columns is not constant")
  if (header)
  {
    nrows <- length(row.lens) - 1
    ncols <- row.lens[2]
    col.names <- scan (file, what = "", sep = sep, nlines = 1, quiet = TRUE, skip = skip)
    x <- scan (file, sep = sep, skip = skip + 1, quiet = T)
  }
  else
  {
    nrows <- length(row.lens)
    ncols <- row.lens[1]
    x <- scan (file, sep = sep, skip = skip, quiet = T)
    col.names <- NULL
  }
  x <- as.double(x)
  if (ncols > 1)
  {
    x <- matrix (x, nc = ncols, nr = nrows, byr = T)
    colnames(x) <- col.names
  }
  else if (ncols == 1)
    x <- as.vector(x)
  else stop ("wrong number of columns")
  return (x)
}

toeplitz <- function (x)
{
  if (!is.vector(x))
    stop ("x is not a vector")
  n <- length (x)
  A <- matrix (0, n, n)
  A <- matrix (x[abs(col(A)-row(A))+1], n, n)
  return (A)
}





