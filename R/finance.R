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
  if (!is.ts(x)) stop ("method is only for time series")
  if (NCOL(x) == 1) stop ("x is not a multivariate time series")
  res <- portfolio.optim.default (as.matrix(x), ...)
  res$px <- ts(res$px,start=start(x),frequency=frequency(x))
  return (res)
}

portfolio.optim.default <- function (x, pm = mean(x), riskless = FALSE, shorts = FALSE, rf = 0.0)
{
  if (!require (quadprog, quietly=TRUE))
    stop ("Package quadprog is needed. Stopping")
  if (NCOL(x) == 1) stop ("x is not a matrix")
  if (any(is.na(x))) stop ("NAs in x")
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

get.hist.quote <- function (instrument = "^gdax", start, end,
                            quote = c("Open", "High", "Low", "Close"), provider = "yahoo",
                            method = c("auto", "wget", "lynx"))
{
  if (!require (chron, quietly=TRUE))
    stop ("Package chron is needed. Stopping")
  quote <- match.arg(quote)
  provider <- match.arg(provider)
  method <- match.arg(method)
  mm <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  if (missing(start)) start <- "1 2 1991"
  if (missing(end)) end <- paste (pmatch(strsplit(date(), " ")[[1]][2],mm),
                                  as.numeric(strsplit(date(), " ")[[1]][3])-1,
                                  strsplit(date(), " ")[[1]][5])
  start <- c(as.numeric(strsplit(start, " ")[[1]][1]),
             as.numeric(strsplit(start, " ")[[1]][2]),
             as.numeric(strsplit(start, " ")[[1]][3]))
  end <- c(as.numeric(strsplit(end, " ")[[1]][1]),
           as.numeric(strsplit(end, " ")[[1]][2]),
           as.numeric(strsplit(end, " ")[[1]][3]))
  if (provider == "yahoo")
  {
    url <- paste ("http://chart.yahoo.com/table.csv?s=", instrument, sep="")
    url <- paste (url, "&a=", start[1], "&b=", start[2], "&c=", start[3], sep="")
    url <- paste (url, "&d=", end[1], "&e=", end[2], "&f=", end[3], sep="")
    url <- paste (url, "&g=d&q=q&y=0&z=", instrument, "&x=.csv", sep="")
    destfile <- tempfile ()
    status <- download.file (url, destfile, method = method)
    if (status != 0) 
    {
      unlink (destfile)
      stop (paste("download error, status", status))
    }
    status <- scan (destfile, "", n=1, sep="\n", quiet=TRUE)
    if (substring(status,1,2) == "No")
    {
      unlink (destfile)
      stop (paste("No data available for", instrument))
    }
    x <- read.table (destfile, header=T, sep=",")
    unlink (destfile)
    nser <- pmatch (quote, c("Open", "High", "Low", "Close")) + 1
    n <- nrow(x)
    ser <- as.vector(x[n:1,nser])
    dat <- dates(gsub("-", " ", as.character(x[n:1,1])), format="day mon y")
    seqdat <- seq.dates (as.numeric(dat)[1],as.numeric(dat)[n])
    idx <- match(as.numeric(dat),as.numeric(seqdat))
    newser <- rep(NA,length(seqdat))
    newser[idx] <- ser
    if (as.numeric(dat)[1] != dates(paste(start[1],start[2],start[3]), format = "m d y"))
      cat (paste("time series starts ", dat[1], "\n", sep=""))
    if (as.numeric(dat)[n] != dates(paste(end[1],end[2],end[3]), format = "m d y"))
      cat (paste("time series ends   ", dat[n], "\n", sep=""))
    return (ts(newser, start = as.numeric(dat)[1], end = as.numeric(dat)[n]))
  }
  else stop ("provider not implemented")
}

