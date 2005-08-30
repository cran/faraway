\name{lawn}
\alias{lawn}
\docType{data}
\title{Cut-off times of lawnmowers}
\description{
Data on the cut-off
times of lawnmowers was collected. 3 machines were randomly selected
from those produced by manufacturers A and B. Each machine was tested
twice at low speed and high speed.
}
\usage{data(lawn)}
\format{
  A data frame with 24 observations on the following 4 variables.
  \describe{
    \item{\code{manufact}}{Manufacturer - a factor with levels \code{A} \code{B}}
    \item{\code{machine}}{Lawn mower - a factor with levels \code{m1} \code{m2} \code{m3} \code{m4} \code{m5} \code{m6}}
    \item{\code{speed}}{Speed of testing - a factor with levels \code{H} \code{L}}
    \item{\code{time}}{cut-off time}
  }
}
\source{
  Unknown.
}
\keyword{datasets}
