\name{bliss}
\alias{bliss}
\docType{data}
\title{
Insect mortality due to insecticide
}
\description{
An experiment measuring death rates for insects, with 30 insects at each of five treatment levels.
}
\usage{data(bliss)}
\format{
  A data frame with 5 observations on the following 3 variables.
  \describe{
    \item{\code{dead}}{number dead}
    \item{\code{alive}}{number alive}
    \item{\code{conc}}{concentration of insecticide}
  }
}
\source{
Bliss (1935). The calculation of the dosage-mortality curve. Annals of Applied Biology 22, 134-167.
}
\examples{
data(bliss)
## maybe str(bliss) ; plot(bliss) ...
}
\keyword{datasets}
