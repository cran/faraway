\name{nels88}
\alias{nels88}
\docType{data}
\title{
National Education Longitudinal Study of 1988 
}
\description{
A subset of the National Education Longitudinal Study of 1988 
}
\usage{data(nels88)}
\format{
  A data frame with 260 observations on the following 5 variables.
  \describe{
    \item{\code{sex}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{race}}{a factor with levels \code{White} \code{Asian} \code{Black} \code{Hispanic}}
    \item{\code{ses}}{a numeric vector}
    \item{\code{paredu}}{a factor with levels \code{ba} \code{college} \code{hs} \code{lesshs} \code{ma} \code{phd}}
    \item{\code{math}}{a numeric vector}
  }
}
\source{
\url{http://www.icpsr.umich.edu/icpsrweb/ICPSR/series/107}
}
\examples{
data(nels88)
## maybe str(nels88) ; plot(nels88) ...
}
\keyword{datasets}
