\name{quadmap}
\alias{quadmap}
\title{Quadratic Map (Logistic Equation)}
\description{
  Computes the quadratic map simulation.
}
\usage{
quadmap(xi = 0.2, a = 4.0, n = 1000)
}
\arguments{
  \item{xi}{the initial value for the iteration.}
  \item{a}{the quadratic map parameter.}
  \item{n}{the length of the simulated series.}
}
\value{
  A vector containing the simulated series.
}
\author{A. Trapletti}
\examples{
x <- quadmap()
acf(x, 10)
}
\keyword{ts}
