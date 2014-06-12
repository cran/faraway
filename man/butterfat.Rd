\name{butterfat}
\alias{butterfat}
\docType{data}
\title{
Butterfat content of milk by breed
}
\description{
Average butterfat content
  (percentages) of milk for random samples of twenty cows (ten
  two-year old and ten mature (greater than four years old)) from each
  of five breeds. The data are from Canadian records of pure-bred
  dairy cattle.
}
\usage{data(butterfat)}
\format{
  A data frame with 100 observations on the following 3 variables.
  \describe{
    \item{\code{Butterfat}}{butter fat content by percentage}
    \item{\code{Breed}}{a factor with levels \code{Ayrshire} \code{Canadian} \code{Guernsey} \code{Holstein-Fresian} \code{Jersey}}
    \item{\code{Age}}{a factor with levels \code{2year} \code{Mature}}
  }
}
\source{
Unknown
}
\examples{
data(butterfat)
## maybe str(butterfat) ; plot(butterfat) ...
}
\keyword{datasets}
