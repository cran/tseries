## Copyright (C) 1997-2001  Adrian Trapletti
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
## Financial time series analysis
##

portfolio.optim <- function (x, ...) UseMethod ("portfolio.optim")

portfolio.optim.ts <-
function (x, ...)
{
    if(!is.ts(x))
        stop("method is only for time series")
    if(NCOL(x) == 1)
        stop("x is not a multivariate time series")
    res <- portfolio.optim.default(as.matrix(x), ...)
    res$px <- ts(res$px, start = start(x), frequency = frequency(x))
    return(res)
}

portfolio.optim.default <-
function(x, pm = mean(x), riskless = FALSE, shorts = FALSE,
         rf = 0.0, reslow = NULL, reshigh = NULL, covmat = cov(x), ...) 
{
    if(!require(quadprog, quietly=TRUE))
        stop("Package quadprog is needed. Stopping")
    if(NCOL(x) == 1)
        stop("x is not a matrix")
    if(any(is.na(x)))
        stop("NAs in x")
    k <- dim(x)[2]
    if(!is.matrix(covmat)) {
        stop("covmat is not a matrix")
    }
    if((dim(covmat)[1] !=k) | (dim(covmat)[2] !=k)) {
      stop("covmat has not the right dimension")
    }
    Dmat <- covmat
    dvec <- rep(0, k)
    big <- 1e+100
    if(!is.null(reslow) & is.null(reshigh)) {
        reshigh <- rep(big, k)
    }
    if(is.null(reslow) & !is.null(reshigh)) {
        reslow <- -rep(big, k)
    }
    if(!is.null(reslow)) {
        if(!is.vector(reslow)) {
            stop("reslow is not a vector")
        }
        if(length(reslow) != k) {
            stop("reslow has not the right dimension")
        }
    }
    if(!is.null(reshigh)) {
        if(!is.vector(reshigh)) {
            stop("reshigh is not a vector")
        }
        if(length(reshigh) != k) {
            stop("reshigh has not the right dimension")
        }
    }
    if(riskless) {
        a1 <- apply(x, 2, mean)-rf
        if(shorts) {
            a2 <- NULL
            b2 <- NULL
        }
        else {
            a2 <- matrix(0, k, k)
            diag(a2) <- 1
            b2 <- rep(0, k)
        }
        if(!is.null(reslow) & !is.null(reshigh)) {
            a3 <- matrix(0, k, k)
            diag(a3) <- 1
            Amat <- t(rbind(a1, a2, a3, -a3))
            b0 <- c(pm-rf, b2, reslow, -reshigh)
        }
        else {
            Amat <- t(rbind(a1, a2))
            b0 <- c(pm-rf, b2)
        }
        res <- solve.QP(Dmat, dvec, Amat, bvec=b0, meq=1)
    }
    else {
        a1 <- rep(1, k)
        a2 <- apply(x, 2, mean)
        if(shorts) {
            if(!is.null(reslow) & !is.null(reshigh)) {
                a3 <- matrix(0, k, k)
                diag(a3) <- 1
                Amat <- t(rbind(a1, a2, a3, -a3))
                b0 <- c(1, pm, reslow, -reshigh)
            }
            else {
                Amat <- t(rbind(a1, a2))
                b0 <- c(1, pm)
            }              
        }
        else {
            a3 <- matrix(0, k, k)
            diag(a3) <- 1
            b3 <- rep(0, k)
            if(!is.null(reslow) & !is.null(reshigh)) {
                Amat <- t(rbind(a1, a2, a3, a3, -a3))
                b0 <- c(1, pm, b3, reslow, -reshigh)
            }
            else {
                Amat <- t(rbind(a1, a2, a3))
                b0 <- c(1, pm, b3)
            }
        }
        res <- solve.QP(Dmat, dvec, Amat, bvec=b0, meq=2)
    }
    y <- t(res$solution%*%t(x))
    ans <- list(pw = res$solution, px = y, pm = mean(y), ps = sd(y))
    return(ans)
}

get.hist.quote <-
function(instrument = "^gdax", start, end,
         quote = c("Open", "High", "Low", "Close") , provider = "yahoo", method = "auto")
{
    if(missing(start)) start <- "1991-01-02"
    if(missing(end)) end <- format(Sys.time() - 86400, "%Y-%m-%d")
  
    provider <- match.arg(provider)

    start <- as.POSIXct(start, tz = "GMT")
    end <- as.POSIXct(end, tz = "GMT")

    if(provider == "yahoo") {
        url <- paste("http://chart.yahoo.com/table.csv?s=",
                     instrument,
                     format(start, "&a=%m&b=%d&c=%Y"),
                     format(end, "&d=%m&e=%d&f=%Y"),
                     "&g=d&q=q&y=0&z=",
                     instrument,
                     "&x=.csv",
                     sep = "")
        destfile <- tempfile()
        status <- download.file(url, destfile, method = method)
        if(status != 0) {
            unlink(destfile)
            stop(paste("download error, status", status))
        }
        status <- scan(destfile, "", n = 1, sep = "\n", quiet = TRUE)
        if(substring(status, 1, 2) == "No") {
            unlink(destfile)
            stop(paste("No data available for", instrument))
        }
        x <- read.table(destfile, header = TRUE, sep = ",")
        unlink(destfile)

        nser <- pmatch(quote, names(x)[-1]) + 1
        if(any(is.na(nser)))
            stop("This quote is not available")
        n <- nrow(x)

        ## Yahoo currently formats dates as `26-Jun-01', hence need C
        ## LC_TIME locale for getting the month right.
        lct <- Sys.getlocale("LC_TIME")
        Sys.setlocale("LC_TIME", "C")
        on.exit(Sys.setlocale("LC_TIME", lct))

        dat <- gsub(" ", "0", as.character(x[, 1])) # Need the gsub?
        dat <- as.POSIXct(strptime(dat, "%d-%b-%y"), tz = "GMT")
        if(dat[n] != start)
            cat(format(dat[n], "time series starts %Y-%m-%d\n"))
        if(dat[1] != end)
            cat(format(dat[1], "time series ends   %Y-%m-%d\n"))
        jdat <- julian(dat)
        ind <- jdat - jdat[n] + 1
        y <- matrix(NA, nr = max(ind), nc = length(nser))
        y[ind, ] <- as.matrix(x[, nser, drop = FALSE])
        colnames(y) <- names(x)[nser]
        return(ts(y, start = jdat[n], end = jdat[1]))
    }
    else stop("provider not implemented")
}

maxdrawdown <-
function(x)
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    cmaxx <- cummax(x)-x
    mdd <- max(cmaxx)
    to <- which(mdd == cmaxx)
    from <- double(NROW(to))
    for (i in 1:NROW(to))
        from[i] <- max(which(cmaxx[1:to[i]] == 0))
    return(list(maxdrawdown = mdd, from = from, to = to))
}

sharpe <-
function(x, r = 0, scale = sqrt(250))
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(NROW(x) == 1)
        return(NA)
    else {
        y <- diff(x)
        return(scale * (mean(y)-r)/sd(y))
    }
}

sterling <-
function(x)
{
    if(NCOL(x) > 1)
        stop("x is not a vector or univariate time series")
    if(any(is.na(x)))
        stop("NAs in x")
    if(NROW(x) == 1)
        return(NA)
    else {
        return((x[NROW(x)]-x[1]) / maxdrawdown(x)$maxdrawdown)
    }
}

plotOHLC <-
function(x, xlim = NULL, ylim = NULL, xlab = "Time", ylab,
         col = par("col"), bg = par("bg"), axes = TRUE,
         frame.plot = axes, ann = par("ann"), main = NULL,
         date = c("calendar", "julian"), format = "%Y-%m-%d", ...)
{
  if ((!is.mts(x)) ||
      (colnames(x)[1] != "Open") ||
      (colnames(x)[2] != "High") ||
      (colnames(x)[3] != "Low") ||
      (colnames(x)[4] != "Close"))
      stop("x is not a open/high/low/close time series")
  xlabel <- if (!missing(x)) 
      deparse(substitute(x))
  else NULL
  if (missing(ylab)) 
      ylab <- xlabel
  date <- match.arg(date)
  time.x <- time(x)
  dt <- min(lag(time.x)-time.x)/3
  ylim <- c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
  if (is.null(xlim)) 
      xlim <- range(time.x)
  if (is.null(ylim)) 
      ylim <- range(x[is.finite(x)])
  plot.new()
  plot.window(xlim, ylim, ...)
  for (i in 1:NROW(x)) {
      segments(time.x[i], x[i,"High"], time.x[i], x[i,"Low"],
               col = col[1], bg = bg)
      segments(time.x[i] - dt, x[i,"Open"], time.x[i], x[i,"Open"],
               col = col[1], bg = bg)
      segments(time.x[i], x[i,"Close"], time.x[i] + dt, x[i,"Close"],
               col = col[1], bg = bg)
  }
  if (ann) 
      title(main = main, xlab = xlab, ylab = ylab, ...)  
  if (axes) {
      if (date == "julian") {
          axis(1, ...)
          axis(2, ...)
      }
      else {
          n <- NROW(x)
          lab.ind <- round(seq(1, n, length=5))
          D <- as.vector(time.x[lab.ind]*86400) + as.POSIXct("1970-01-01", tz = "GMT")
          DD <- format.POSIXct(D, format = format, tz ="GMT")
          axis(1, at=time.x[lab.ind], lab=DD, ...)
          axis(2, ...)
      }
  }
  if (frame.plot) 
      box(...)
}
