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
## Financial time series analysis
##


portfolio.optim <- function (x, ...) UseMethod ("portfolio.optim")

portfolio.optim.ts <- function (x, ...)
{
    if(!is.ts(x)) stop("method is only for time series")
    if(NCOL(x) == 1) stop("x is not a multivariate time series")
    res <- portfolio.optim.default(as.matrix(x), ...)
    res$px <- ts(res$px, start = start(x), frequency = frequency(x))
    return(res)
}

portfolio.optim.default <-
function(x, pm = mean(x), riskless = FALSE, shorts = FALSE, rf = 0.0)
{
    if(!require(quadprog, quietly=TRUE))
        stop("Package quadprog is needed. Stopping")
    if(NCOL(x) == 1) stop("x is not a matrix")
    if(any(is.na(x))) stop("NAs in x")
    k <- dim(x)[2]
    Dmat <- cov(x)
    dvec <- rep(0, k)
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
        Amat <- t(rbind(a1, a2))
        b0 <- c(pm-rf, b2)
        res <- solve.QP(Dmat, dvec, Amat, bvec=b0, meq=1)
    }
    else {
        a1 <- rep(1, k)
        a2 <- apply(x, 2, mean)
        if(shorts) {
            a3 <- NULL
            b3 <- NULL
        }
        else {
            a3 <- matrix(0, k, k)
            diag(a3) <- 1
            b3 <- rep(0, k)
        }
        Amat <- t(rbind(a1, a2, a3))
        b0 <- c(1, pm, b3)
        res <- solve.QP(Dmat, dvec, Amat, bvec=b0, meq=2)
    }
    y <- t(res$solution%*%t(x))
    ans <- list(pw = res$solution, px = y, pm = mean(y), ps = sd(y))
    return(ans)
}

get.hist.quote <-
function(instrument = "^gdax", start, end,
         quote = "Open", provider = "yahoo", method = "auto")
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
