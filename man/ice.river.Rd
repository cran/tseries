\name{ice.river}
\alias{ice.river}
\alias{flow.vat}
\alias{flow.jok}
\alias{prec}
\alias{temp}
\title{
  Icelandic River Data
}
\description{
  Contains the Icelandic river data as presented in
  \bibcitet{|Tong:1990|pages 432--440}.
}
\usage{
data(ice.river)
}
\format{
  4 univariate time series \code{flow.vat}, \code{flow.jok},
  \code{prec}, and \code{temp}, each with 1095 observations and the
  joint series \code{ice.river}.
} 
\details{
  The series are daily observations from Jan. 1, 1972 to Dec. 31, 1974
  on 4 variables: \code{flow.vat}, mean daily flow of Vatnsdalsa river
  (cms), \code{flow.jok}, mean daily flow of Jokulsa Eystri river (cms),
  \code{prec}, daily precipitation in Hveravellir (mm), and mean daily
  temperature in Hveravellir (deg C).

  These datasets were introduced into the literature in a paper by
  \bibcitet{Tong+Thanoon+Gudmundsson:1985}.
}
\source{
  Time Series Data Library:
  \url{https://robjhyndman.com/TSDL/}
}
\references{
  \bibshow{Tong:1990, Tong+Thanoon+Gudmundsson:1985}
}
\keyword{datasets}
