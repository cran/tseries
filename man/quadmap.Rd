\name{quadmap}
\title{Quadratic Map (Logistic Equation)}
\usage{
quadmap (xi = 0.2, a = 4.0, n = 1000)
}
\alias{quadmap}
\arguments{
  \item{xi}{the initial value for the iteration.}
  \item{a}{the quadratic map parameter.}
  \item{n}{the length of the simulated series.}
}
\description{
  Computes the quadratic map simulation.
}
\value{
  A vector containing the simulated series.
}
\author{A. Trapletti}
\examples{
x <- quadmap ()
acf (x, 10)
}
\keyword{ts}

