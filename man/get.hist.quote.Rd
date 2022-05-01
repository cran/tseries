\name{get.hist.quote}
\alias{get.hist.quote}
\title{Download Historical Finance Data}
\description{
  Download historical financial data from a given data provider over
  the WWW.
}
\usage{
get.hist.quote(instrument = "^gdax", start, end,
               quote = c("Open", "High", "Low", "Close"),
               provider = c("yahoo"), method = NULL,
               origin = "1899-12-30", compression = "d",
	       retclass = c("zoo", "ts"), quiet = FALSE, drop = FALSE) 
}
\arguments{
  \item{instrument}{a character string giving the name of the quote
    symbol to download.  See the web page of the data provider for
    information about the available quote symbols.}
  \item{start}{an R object specifying the date of the start of the
    period to download.  This must be in a form which is recognized by
    \code{\link{as.POSIXct}}, which includes R POSIX date/time objects,
    objects of class \code{"date"} (from package \code{date}) and
    \code{"chron"} and \code{"dates"} (from package \code{chron}), and
    character strings representing dates in ISO 8601 format.  Defaults
    to 1992-01-02.}
  \item{end}{an R object specifying the end of the download period, see
    above.  Defaults to yesterday.}
  \item{quote}{a character string or vector indicating whether to
    download opening, high, low, or closing quotes, or volume.  For the
    default provider, this can be specified as \code{"Open"},
    \code{"High"}, \code{"Low"}, \code{"Close"}, \code{"Adjusted"}, and
    \code{"Volume"}, respectively.  Abbreviations are allowed.}
  \item{provider}{a character string with the name of the data
    provider.  Currently, only \code{"yahoo"} is supported via 
    \code{\link[quantmod:getSymbols.yahoo]{getSymbols}} from package
    \pkg{quantmod} for the Yahoo Finance source.
    Provider \code{"oanda"} is no longer available.}
  \item{method}{No longer used.}
  \item{origin}{an R object specifying the origin of the Julian dates, see
    above.  Defaults to 1899-12-30 (Popular spreadsheet programs
    internally also use Julian dates with this origin).}
  \item{compression}{Governs the granularity of the retrieved data;
    \code{"d"} for daily, \code{"w"} for weekly or \code{"m"} for
    monthly.  Defaults to \code{"d"}.  For the provider \code{"oanda"},
    this argument is ignored.}
  \item{retclass}{character specifying which class the return value
    should have: can be either \code{"zoo"} (with \code{"Date"} index),
    or \code{"ts"} (with numeric index corresponding to days since
    \code{origin}).}
  \item{quiet}{logical. Should status messages (if any) be suppressed?}
  \item{drop}{logical.  If \code{TRUE} the result is coerced to the
    lowest possible dimension.  Default is \code{FALSE}.} 
}
\value{
  A time series containing the data either as a \code{"zoo"} series
  (default) or a \code{"ts"} series.  The \code{"zoo"} series is created
  with \code{\link[zoo]{zoo}} and has an index of class \code{"Date"}.
  If a \code{"ts"} series is returned, the index is in physical time,
  i.e., weekends, holidays, and missing days are filled with \code{NA}s
  if not available.  The time scale is given in Julian dates (days since
  the \code{origin}).
}
\author{A. Trapletti}
\seealso{
  \code{\link[quantmod]{getSymbols}} for downloads from various
  providers; 
  \code{\link[zoo]{zoo}},
  \code{\link{ts}},
  \code{\link{as.Date}},
  \code{\link{as.POSIXct}},
}
\examples{
tryCatch({
  x <- get.hist.quote(instrument = "^gspc", start = "1998-01-01",
                      quote = "Close")
  plot(x)

  x <- get.hist.quote(instrument = "ibm", quote = c("Cl", "Vol"))
  plot(x, main = "International Business Machines Corp")

  spc <- get.hist.quote(instrument = "^gspc", start = "1998-01-01",
         quote = "Close")
  ibm <- get.hist.quote(instrument = "ibm",  start = "1998-01-01",
         quote = "Adj")
  require("zoo")		# For merge() method.
  x <- merge(spc, ibm)
  plot(x, main = "IBM vs S&P 500")
}, error = identity)
}
\keyword{ts}
