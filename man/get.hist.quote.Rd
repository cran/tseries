\name{get.hist.quote}
\title{Download Historical Finance Data}
\usage{
get.hist.quote (instrument = "^gdax", start, end,
                quote = c("Open", "High", "Low", "Close", "Volume"),
                provider = "yahoo",
                method = "auto")
}
\alias{get.hist.quote}
\arguments{
  \item{instrument}{the name of the quote symbol to download. See the
    web page of the data provider for information about the quote symbols.}
  \item{start, end}{strings representing the start and end date
    of the period to download. The format is \code{"mm dd yyyy"}.
    Default values are \code{"1 2 1991"} and yesterday.}
  \item{quote}{download opening-, high-, low-, closing-quotes, or volume?}
  \item{provider}{the name of the data provider. Currently implemented
    is \code{"yahoo"}. See \url{http://quote.yahoo.com/} for more
    information.}
  \item{method}{tool to be used for downloading the data. See
    \code{\link{download.file}} for the available download methods.}
}
\description{
  Download historical financial data from a given data provider over
  the WWW.
}
\value{
  A time series containing the data in physical time, i.e., weekends,
  holidays, and missing days are filled with \code{NA}s. The time scale
  is given in Julian dates.
}
\author{A. Trapletti}
\seealso{
  \code{\link{ts}}, \code{\link{chron}}, \code{\link{download.file}},
  \url{http://quote.yahoo.com/} 
}
\examples{
\dontrun{x <- get.hist.quote (instrument = "^spc", start = "1 1 1998",
                     quote = "Close")
plot (x)}
}
\keyword{ts}

