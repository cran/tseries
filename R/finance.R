## Copyright (C) 1997-2003  Adrian Trapletti
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
    if(!require("quadprog", quietly=TRUE))
        stop("package", sQuote("quadprog"), "is needed.  Stopping")
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
function (instrument = "^gdax", start, end,
          quote = c("Open", "High", "Low", "Close"),
          provider = c("yahoo", "oanda"), method = "auto",
          origin = "1899-12-30", compression = "d",
	  retclass = c("zoo", "its", "ts"))
    ## Added new argument 'compression'.
    ## May be "d", "w" or "m", for daily weekly or monthly.
    ## Defaults to "d".
    ## John Bollinger, 2004-10-27, www.BollingerBands.com, bbands@yahoo.com
    ##
    ## Changed POSIXct class to Date class, 2005-03-31
{
    if(missing(start)) start <- "1991-01-02"
    if(missing(end)) end <- format(Sys.Date() - 1, "%Y-%m-%d")
  
    provider <- match.arg(provider)
    retclass <- match.arg(retclass)

    start <- as.Date(start)
    end <- as.Date(end)

    if(provider == "yahoo") {
        url <-
            paste("http://chart.yahoo.com/table.csv?s=",
                  instrument, 
                  format(start,
                         paste("&a=",
                               as.character(as.numeric(format(start, "%m"))-1),
                               "&b=%d&c=%Y",
                               sep = "")),
                  format(end,
                         paste("&d=",
                               as.character(as.numeric(format(end, "%m"))-1),
                               "&e=%d&f=%Y",
                               sep = "")), 
                  "&g=", compression,
                  "&q=q&y=0&z=", instrument,
                  "&x=.csv",
                  sep = "")
        destfile <- tempfile()
        status <- download.file(url, destfile, method = method)
        if(status != 0) {
            unlink(destfile)
            stop(paste("download error, status", status))
        }
        nlines <- length(count.fields(destfile, sep = "\n"))
        if(nlines == 1) {
            unlink(destfile)
            stop(paste("no data available for", instrument))
        }
        
        ## Yahoo includes rows concerning dividends,
        ## hence need fill = TRUE and na.omit
        x <- read.table(destfile, header = TRUE, sep = ",", as.is = TRUE, fill = TRUE)
        x <- na.omit(x)
        
        ## Debug
        ## cat("read.table: start =", x[NROW(x),"Date"], "\n")
        ## cat("read.table: end   =", x[1,"Date"], "\n")
        
        unlink(destfile)
        
        names(x) <- gsub("\\.", "", names(x))
        nser <- pmatch(quote, names(x)[-1]) + 1
        if(any(is.na(nser)))
            stop("this quote is not available")
        n <- nrow(x)

        ## Yahoo currently formats dates as '26-Jun-01', hence need C
        ## LC_TIME locale for getting the month right.
        lct <- Sys.getlocale("LC_TIME")
        Sys.setlocale("LC_TIME", "C")
        on.exit(Sys.setlocale("LC_TIME", lct))

        dat <- gsub(" ", "0", as.character(x[, 1])) # Need the gsub?
        idx <- c(grep(".*-0.", dat),  grep(".*-1.", dat))
        dat[idx] <- paste(substr(dat[idx], 1, nchar(dat[idx]) - 2),
                          "20",
                          substr(dat[idx], nchar(dat[idx]) - 1, nchar(dat[idx])),
                          sep = "")
        dat[-idx] <- paste(substr(dat[-idx], 1, nchar(dat[-idx]) - 2),
                           "19",
                           substr(dat[-idx], nchar(dat[-idx]) - 1, nchar(dat[-idx])),
                           sep = "")
        dat <- as.Date(dat, "%d-%b-%Y")
        if(dat[n] != start)
            cat(format(dat[n], "time series starts %Y-%m-%d\n"))
        if(dat[1] != end)
            cat(format(dat[1], "time series ends   %Y-%m-%d\n"))

	if(retclass == "ts") {
            jdat <- unclass(julian(dat, origin = as.Date(origin)))
            ## We need unclass() because 1.7.0 does not allow adding a number
            ## to a "difftime" object. 
            ind <- jdat - jdat[n] + 1
            y <- matrix(NA, nr = max(ind), nc = length(nser))
            y[ind, ] <- as.matrix(x[, nser, drop = FALSE])
            colnames(y) <- names(x)[nser]
            return(ts(y, start = jdat[n], end = jdat[1]))
	} else {
	  x <- as.matrix(x[, nser, drop = FALSE])
	  rownames(x) <- NULL
	  y <- zoo(x, dat)
	  if(retclass == "its") {
	    if("package:its" %in% search() || require("its", quietly = TRUE)) {
	        index(y) <- as.POSIXct(index(y))
	        y <- as.its(y)
	    } else {
	      warning("package its could not be loaded: zoo series returned")
	    }
	  }
	  return(y)
	}
    }
    else if(provider == "oanda") {
        if(!missing(quote)) {
            warning("argument 'quote' ignored for provider 'oanda'")
        }
        if(!missing(compression)) {
            warning("argument 'compression' ignored for provider 'oanda'")
        }
        
        url <-
            paste("http://www.oanda.com/convert/fxhistory?lang=en&date1=",
                  format(start, "%m"), "%2F", format(start, "%d"), "%2F", format(start, "%Y"),
                  "&date=",
                  format(end, "%m"), "%2F", format(end, "%d"), "%2F", format(end, "%Y"),
                  "&date_fmt=us&exch=",
                  unlist(strsplit(instrument, split = "/"))[1],
                  "&exch2=&expr=",
                  unlist(strsplit(instrument, split = "/"))[2],
                  "&expr2=&margin_fixed=0&&SUBMIT=Get+Table&format=ASCII&redirected=1",
                  sep="")
        destfile <- tempfile()
        
        status <- download.file(url, destfile, method = method)
        if(status != 0) {
            unlink(destfile)
            stop(paste("download error, status", status))
        }
        
        x <- readLines(destfile)
        unlink(destfile)
        
        if(length(grep("Sorry", x)) > 0) {
            msg <- unlist(strsplit(gsub("<[a-zA-Z0-9\\/]*>", "", x[grep("Sorry", x)]), split = " "))
            msg <- paste(msg[msg != ""], collapse = " ")
            stop("Message from Oanda: ", msg)
        }
        
        first <- which(substr(x, 1, 5) == "<PRE>")
        last <- which(x == "</PRE>") - 1
        if((length(first) == 0) || (length(last) == 0)) {
            stop(paste("no data available for", instrument))
        }

        x[first] <- substr(x[first], 6, nchar(x[first]))
        split <- strsplit(x[first:last], split = " ")
        x <- cbind(unlist(lapply(split, function(x) { x[1] })), unlist(lapply(split, function(x) { x[length(x)] })))
        n <- nrow(x)
        
        dat <- as.Date(x[,1], format = "%m/%d/%Y")
        if(dat[1] != start)
            cat(format(dat[1], "time series starts %Y-%m-%d\n"))
        if(dat[n] != end)
            cat(format(dat[n], "time series ends   %Y-%m-%d\n"))

	if(retclass == "ts") {
            jdat <- unclass(julian(dat, origin = as.Date(origin)))       
            ind <- jdat - jdat[1] + 1
            y <- rep(NA, max(ind))
            y[ind] <- as.numeric(x[,2])
            return(ts(y, start = jdat[1], end = jdat[n]))
	} else {	  
	  y <- zoo(as.numeric(x[,2]), dat)
	  if(retclass == "its") {
	    if("package:its" %in% search() || require("its", quietly = TRUE)) {
	        index(y) <- as.POSIXct(index(y))
	        y <- as.its(y)
	    } else {
	      warning("package its could not be loaded: zoo series returned")
	    }
	  }
	  return(y)
	}
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
         date = c("calendar", "julian"), format = "%Y-%m-%d",
         origin = "1899-12-30", ...)
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
  if (is.null(xlim)) 
      xlim <- range(time.x)
  if (is.null(ylim)) 
      ylim <- range(x[is.finite(x)])
  plot.new()
  plot.window(xlim, ylim, ...)
  segments(time.x, x[, "High"], time.x, x[, "Low"], col = col[1], bg = bg)
  segments(time.x - dt, x[, "Open"], time.x, x[, "Open"],
           col = col[1], bg = bg)
  segments(time.x, x[, "Close"], time.x + dt, x[, "Close"],
           col = col[1], bg = bg)
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
          D <- as.vector(time.x[lab.ind]*86400) + as.POSIXct(origin, tz = "GMT")
          DD <- format.POSIXct(D, format = format, tz ="GMT")
          axis(1, at=time.x[lab.ind], lab=DD, ...)
          axis(2, ...)
      }
  }
  if (frame.plot) 
      box(...)
}

