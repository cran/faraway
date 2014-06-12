\name{punting}
\alias{punting}
\docType{data}
\title{
Leg strength and punting
}
\description{
 Investigators studied physical characteristics and ability in 13
  (American) football punters. Each volunteer punted a football ten times. The
  investigators recorded the average distance for the ten punts, in
  feet.
}
\usage{data(punting)}
\format{
  A data frame with 13 observations on the following 7 variables.
  \describe{
    \item{\code{Distance}}{average distance over 10 punts}
    \item{\code{Hang}}{hang time}
    \item{\code{RStr}}{right leg strength in pounds}
    \item{\code{LStr}}{left leg strength in pounds}
    \item{\code{RFlex}}{right hamstring muscle flexibility in degrees}
    \item{\code{LFlex}}{left hamstring muscle flexibility in degrees}
    \item{\code{OStr}}{overall leg strength in foot pounds}
  }
}
\source{
Unknown
}
\examples{
data(punting)
## maybe str(punting) ; plot(punting) ...
}
\keyword{datasets}
