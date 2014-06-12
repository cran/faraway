\name{insulgas}
\alias{insulgas}
\docType{data}
\title{
Effects of insulation on gas consumption
}
\description{
Data on natural gas usage in a house.
  The weekly gas consumption (in 1000 cubic feet) and the average
  outside temperature (in degrees Celsius) was recorded for 26 weeks
  before and 30 weeks after cavity-wall insulation had been installed.
  The house thermostat was set at 20C throughout.
}
\usage{data(insulgas)}
\format{
  A data frame with 44 observations on the following 3 variables.
  \describe{
    \item{\code{Insulate}}{a factor with levels \code{After} \code{Before}}
    \item{\code{Temp}}{Outside temperature}
    \item{\code{Gas}}{Weekly consumption in 1000 cubic feet}
  }
}
\source{
MASS package as whiteside
}
\examples{
data(insulgas)
## maybe str(insulgas) ; plot(insulgas) ...
}
\keyword{datasets}
