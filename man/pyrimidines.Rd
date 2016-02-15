\name{pyrimidines}
\alias{pyrimidines}
\docType{data}
\title{
Activity in pyrimidines
}
\description{
  Structural information on 74 2,4-diamino- 5-(substituted benzyl)
  pyrimidines used as inhibitors of DHFR in E. coli. There are 3
  positions where chemical activity occurs and 9 attributes per position
  leading to 27 total predictors. One predictor had no variability and was removed from the data set.
26 chemical properties of 74 compounds and an activity level
}
\usage{data(pyrimidines)}
\format{
  A data frame with 74 observations on the following 27 variables.
  \describe{
    \item{\code{p1.polar}}{measured on a [0,1] scale}
    \item{\code{p1.size}}{measured on a [0,1] scale}
    \item{\code{p1.flex}}{measured on a [0,1] scale}
    \item{\code{p1.h.doner}}{measured on a [0,1] scale}
    \item{\code{p1.h.acceptor}}{measured on a [0,1] scale}
    \item{\code{p1.pi.doner}}{measured on a [0,1] scale}
    \item{\code{p1.pi.acceptor}}{measured on a [0,1] scale}
    \item{\code{p1.polarisable}}{measured on a [0,1] scale}
    \item{\code{p1.sigma}}{measured on a [0,1] scale}
    \item{\code{p2.polar}}{measured on a [0,1] scale}
    \item{\code{p2.size}}{measured on a [0,1] scale}
    \item{\code{p2.flex}}{measured on a [0,1] scale}
    \item{\code{p2.h.doner}}{measured on a [0,1] scale}
    \item{\code{p2.h.acceptor}}{measured on a [0,1] scale}
    \item{\code{p2.pi.doner}}{measured on a [0,1] scale}
    \item{\code{p2.pi.acceptor}}{measured on a [0,1] scale}
    \item{\code{p2.polarisable}}{measured on a [0,1] scale}
    \item{\code{p2.sigma}}{measured on a [0,1] scale}
    \item{\code{p3.polar}}{measured on a [0,1] scale}
    \item{\code{p3.size}}{measured on a [0,1] scale}
    \item{\code{p3.flex}}{measured on a [0,1] scale}
    \item{\code{p3.h.doner}}{measured on a [0,1] scale}
    \item{\code{p3.h.acceptor}}{measured on a [0,1] scale}
    \item{\code{p3.pi.doner}}{measured on a [0,1] scale}
    \item{\code{p3.polarisable}}{measured on a [0,1] scale}
    \item{\code{p3.sigma}}{measured on a [0,1] scale}
    \item{\code{activity}}{log 1/Ki, where Ki is the inhibition constant
      as experimentally assayed, scaled to [0,1]}
  }
}
\source{
 Jonathan D. Hirst, Ross D. King, Michael J. E. Sternberg (1994) Quantitative structure-activity relationships by neural networks and inductive logic programming. I. The inhibition of dihydrofolate reductase by pyrimidines
 \url{http://dx.doi.org/10.1007/BF00125375}
 }
\examples{
data(pyrimidines)
## maybe str(pyrimidines) ; plot(pyrimidines) ...
}
\keyword{datasets}
