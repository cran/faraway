\name{penicillin}
\alias{penicillin}
\docType{data}
\title{
Penicillin yield by block and treatment
}
\description{
The production of
  penicillin uses a raw material, corn steep liquor, is quite
  variable and can only be made in blends sufficient for four
  runs. There are four processes, A, B, C and D, for the production.
}
\usage{data(penicillin)}
\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{\code{treat}}{a factor with levels \code{A} \code{B} \code{C} \code{D}}
    \item{\code{blend}}{a factor with levels \code{Blend1} \code{Blend2} \code{Blend3} \code{Blend4} \code{Blend5}}
    \item{\code{yield}}{a numeric vector}
  }
}
\source{
Box, G., W. Hunter, and J. Hunter (1978). Statistics for Experimenters. New York: Wiley.
}
\examples{
data(penicillin)
## maybe str(penicillin) ; plot(penicillin) ...
}
\keyword{datasets}
