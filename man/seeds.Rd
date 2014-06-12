\name{seeds}
\alias{seeds}
\docType{data}
\title{
Germination of seeds depending on moisture and covering
}
\description{
A Biologist analyzed an experiment to
determine the effect of moisture content on seed
germination. Eight boxes of 100 seeds each were
treated with the same moisture level. 4 boxes were
covered and 4 left uncovered. The process was repeated
at 6 different moisture levels (nonlinear scale).
}
\usage{data(seeds)}
\format{
  A data frame with 48 observations on the following 3 variables.
  \describe{
    \item{\code{germ}}{percentage germinated}
    \item{\code{moisture}}{moisture level}
    \item{\code{covered}}{a factor with levels \code{no} \code{yes}}
  }
}
\source{
Chatfield C. (1995) Problem Solving: A Statistician's Guide, 2ed Chapman Hall.
}
\examples{
data(seeds)
## maybe str(seeds) ; plot(seeds) ...
}
\keyword{datasets}
