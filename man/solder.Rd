\name{solder}
\alias{solder}
\docType{data}
\title{
Solder skips in printing circuit boards
}
\description{
ATT ran an experiment varying five factors
relevant to a wave-soldering procedure for mounting components on
printed circuit boards.  The response variable, skips, is a count of
how many solder skips appeared to a visual inspection.
}
\usage{data(solder)}
\format{
  A data frame with 900 observations on the following 6 variables.
  \describe{
    \item{\code{Opening}}{a factor with levels \code{L} \code{M} \code{S}}
    \item{\code{Solder}}{a factor with levels \code{Thick} \code{Thin}}
    \item{\code{Mask}}{a factor with levels \code{A1.5} \code{A3} \code{A6} \code{B3} \code{B6}}
    \item{\code{PadType}}{a factor with levels \code{D4} \code{D6} \code{D7} \code{L4} \code{L6} \code{L7} \code{L8} \code{L9} \code{W4} \code{W9}}
    \item{\code{Panel}}{a numeric vector}
    \item{\code{skips}}{count of how many solder skips appeared to a visual inspection}
  }
}
\source{
Comizzoli, R. B., J. M. Landwehr, and J. D. Sinclair (1990). Robust materials and
processes: Key to reliability. AT&T Technical Journal 69(6), 113-128.
}
\examples{
data(solder)
## maybe str(solder) ; plot(solder) ...
}
\keyword{datasets}
