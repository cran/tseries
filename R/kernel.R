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
# Smooting kernel object
#


kernel <- function (coef, m, name)
{
  if (!is.vector(coef))
    stop ("coef must be a vector")
  if ((length(coef) != m+1) | (length(coef) <= 0))
    stop ("coef has not the correct length")
  kernel <- list (coef=coef, m=m)
  attr(kernel, "name") <- name
  class(kernel) <- "kernel"
  eps <- 1.0e-10
  if ((sum(kernel[-m:m]) > 1.0+eps) || (sum(kernel[-m:m]) < 1.0-eps))
    stop ("coefficients do not add to 1")
  return (kernel) 
}

print.kernel <- function (k, digits = max(3,.Options$digits-3))
{
  cat (attr(k,"name"),"\n")
  for (i in (-k$m:k$m))
    cat (c("coef[",formatC(i,format="d",flag=" "),"] = ",
           formatC(k[i],digits=digits,format="f")),"\n")
}

plot.kernel <- function (k)
{
  y <- c(rev(k$coef[2:(k$m+1)]),k$coef)
  plot ((-k$m:k$m),y,xlab="k",ylab="W[k]",type="h",main=attr(k,"name"))
}

daniell.kernel <- function (m)
{
  if (m < 0) stop ("m is negative")
  return (kernel(rep(1/(2*m+1),m+1),m,paste("Daniell(",m,")",sep="")))
}

fejer.kernel <- function (m, r)
{
  if (r < 1)
    stop ("r is less than 1")
  if (m < 1)
    stop ("m is less than 1")
  n <- 2*m+1
  wn <- double(m+1)
  for (j in (1:m))
  {
    wj <- 2*pi*j/n
    wn[j+1] <- sin(r*wj/2)^2/sin(wj/2)^2
  }
  wn <- wn/(n*r)
  wn[1] <- r/n
  wn <- wn/sum(c(rev(wn[2:(m+1)]),wn))
  return (kernel(wn,m,paste("Fejer(",m,",",r,")",sep="")))
}

df.kernel <- function (k)
{ 
  return (2/sum(k[-k$m:k$m]^2))
}

dirichlet.kernel <- function (m, r)
{
  if (r < 0)
    stop ("r is less than 0")
  if (m < 1)
    stop ("m is less than 1")
  n <- 2*m+1
  wn <- double(m+1)
  for (j in (1:m))
  {
    wj <- 2*pi*j/n
    wn[j+1] <- sin((r+0.5)*wj)/sin(wj/2)
  }
  wn <- wn/n
  wn[1] <- (2*r+1)/n
  wn <- wn/sum(c(rev(wn[2:(m+1)]),wn))
  return (kernel(wn,m,paste("Dirichlet(",m,",",r,")",sep="")))
}

"[.kernel" <- function (k, i)
{
  y <- c(rev(k$coef[2:(k$m+1)]),k$coef)
  return (y[i+k$m+1])
}

is.kernel <- function (k)
{
  return (inherits(k, "kernel"))
}

apply.kernel <- function (x, ...)
{
  UseMethod("apply.kernel")
}

apply.kernel.vector <- function (x, k, circular = FALSE)
{
  if (!is.vector(x))
    stop ("x is not a vector")
  if (!is.kernel(k))
    stop ("k is not a kernel")
  if (length(x) <= 2*k$m)
    stop ("x is shorter than kernel k") 
  if (k$m == 0)
    return (x)
  else
  {
    n <- length(x)
    w <- c(k[0:k$m],rep(0,n-2*k$m-1),k[-k$m:-1])
    y <- fft(fft(x)*fft(w),inverse=T)/n
    if (is.numeric(x)) y <- Re(y)
    if (circular)
      return (y)
    else
      return (y[(1+k$m):(n-k$m)])
  }
}

apply.kernel.default <- function (x, k, circular = FALSE)
{  
  if (is.vector(x))
    return (apply.kernel.vector(x,k,circular=circular))
  else if (is.matrix(x))
    return (apply(x,MARGIN=2,FUN=apply.kernel,k,circular=circular))
  else
    stop ("apply.kernel is not available for object x")
}

apply.kernel.ts <- function (x, k, circular = FALSE)
{  
  if (is.univariate.ts(x))
    y <- apply.kernel.vector(as.vector(x),k,circular=circular)
  else 
    y <- apply(x,MARGIN=2,FUN=apply.kernel,k,circular=circular)
  return (ts (y, end=end(x), frequency=frequency(x)))
}

apply.kernel.kernel <- function (k1, k2)
{
  if (!is.kernel(k1))
    stop ("k1 is not a kernel")
  if (!is.kernel(k2))
    stop ("k2 is not a kernel")
  n <- k2$m
  x <- c(rep(0,n),k1[-k1$m:k1$m],rep(0,n))
  coef <- apply.kernel(x,k2,circular=T)
  m <- length(coef)%/%2
  return (kernel(coef[(m+1):length(coef)],m,
                 paste("Composite(", attr(k1, "name"),",",attr(k2, "name"),")",sep="")))
}







