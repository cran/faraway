\name{fround}
\alias{fround}
\alias{pfround}

\title{Formating the Rounding of Numbers}

\description{
    \code{fround} rounds the values in its first argument to the specified
    number of decimal places with surrounding quotes.

    \code{pfround} rounds the values in its first argument to the specified
    number of decimal places without surrounding quotes.

}

\usage{
fround(x, digits)
pfround(x, digits)
}

\arguments{
  \item{x}{a numeric vector.}
  \item{digits}{integer indicating the precision to be used.}
}

\references{Copied from the \code{arm} package}


\author{Andrew Gelman; Yu-Sung Su}

\seealso{
    \code{\link{round}}
        }

\examples{
    x <- 3.1415926
    fround(x, digits=2)
    pfround(x, digits=2)
}

\keyword{manip}
\keyword{print}
