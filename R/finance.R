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
# Financial time series analysis
#


portfolio.optim <- function (obj, ...) { UseMethod ("portfolio.optim") }

portfolio.optim.ts <- function (x, ...)
{
  if (is.univariate.ts(x)) stop ("x is a univariate time series")
  res <- portfolio.optim.default (as.matrix(x), ...)
  res$px <- ts(res$px,start=start(x),frequency=frequency(x))
  return (res)
}

portfolio.optim.default <- function (x, pm = mean(x), riskless = FALSE, shorts = FALSE, rf = 0.0)
{
  if (!require (quadprog, quietly=T)) stop ("Stopping")
  if (!is.matrix(x)) stop ("x is not a matrix")
  k <- dim(x)[2]
  Dmat <- cov(x)
  dvec <- rep(0,k)
  if (riskless)
  {
    a1 <- apply(x,2,mean)-rf
    if (shorts)
    {
      a2 <- NULL
      b2 <- NULL
    }
    else
    {
      a2 <- matrix(0,k,k)
      diag(a2) <- 1
      b2 <- rep(0,k)
    }
    Amat <- t(rbind(a1,a2))
    b0 <- c(pm-rf,b2)
    res <- solve.QP(Dmat,dvec,Amat,bvec=b0,meq=1)
  }
  else
  {
    a1 <- rep(1,k)
    a2 <- apply(x,2,mean)
    if (shorts)
    {
      a3 <- NULL
      b3 <- NULL
    }
    else
    {
      a3 <- matrix(0,k,k)
      diag(a3) <- 1
      b3 <- rep(0,k)
    }  
    Amat <- t(rbind(a1,a2,a3))
    b0 <- c(1,pm,b3)
    res <- solve.QP(Dmat,dvec,Amat,bvec=b0,meq=2)
  }
  y <- t(res$solution%*%t(x))
  ans <- list (pw=res$solution, px=y, pm=mean(y), ps=sd(y))
  return (ans)
}
