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
               provider = "yahoo", method = "auto",
               origin = "1899-12-30", compression = "d") 
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
    \code{"High"}, \code{"Low"}, \code{"Close"}, and \code{"Volume"},
    respectively.  Abbreviations are allowed.}
  \item{provider}{a character string with the name of the data
    provider.  Currently, only \code{"yahoo"} is implemented.  See
    \url{http://quote.yahoo.com/} for more information.}
  \item{method}{tool to be used for downloading the data.  See
    \code{\link{download.file}} for the available download methods.}
  \item{origin}{an R object specifying the origin of the Julian dates, see
    above.  Defaults to 1899-12-30 (Popular spreadsheet programs
    internally also use Julian dates with this origin).}
  \item{compression}{Governs the granularity of the retrieved data;
    \code{"d"} for daily, \code{"w"} for weekly or \code{"m"} for
    monthly.  Defaults to \code{"d"}.}
}
\value{
  A time series containing the data in physical time, i.e., weekends,
  holidays, and missing days are filled with \code{NA}s.  The time scale
  is given in Julian dates (days since the \code{origin}). 
}
\author{A. Trapletti}
\seealso{
  \code{\link{ts}},
  \code{\link{as.POSIXct}},
  \code{\link{download.file}};
  \url{http://quote.yahoo.com/} 
}
\examples{
if(!inherits(try(open(url("http://quote.yahoo.com")), silent = TRUE),
             "try-error")) {
  x <- get.hist.quote(instrument = "^gspc", start = "1998-01-01",
                      quote = "Close")
  plot(x)

  x <- get.hist.quote(instrument = "ibm", quote = c("Cl", "Vol"))
  plot(x, main = "International Business Machines Corp")

  spc <- get.hist.quote(instrument = "^gspc", start = "1998-01-01")
  ibm <- get.hist.quote(instrument = "ibm",  start = "1998-01-01")
  x <- na.remove(ts.union(spc, ibm))
  plot(x, main = "IBM vs S&P 500")
}
}
\keyword{ts}
