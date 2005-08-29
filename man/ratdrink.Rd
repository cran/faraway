\name{ratdrink}
\alias{ratdrink}
\docType{data}
\title{Rat growth weights affected by additives}
\description{
The data
consist of 5 weekly measurements of body weight for 27 rats. The
first 10 rats are on a control treatment while 7 rats have thyroxine
added to their drinking water. 10 Rats have thiouracil added to their water.
}
\usage{data(ratdrink)}
\format{
  A data frame with 135 observations on the following 4 variables.
  \describe{
    \item{\code{wt}}{Weight of the rat}
    \item{\code{weeks}}{Week of the study from 0 to 4}
    \item{\code{subject}}{the rat code number}
    \item{\code{treat}}{treatment applied to the rat drinking water - a factor with levels \code{control} \code{thiouracil} \code{thyroxine}}
  }
}
\source{
  Unknown
}
\keyword{datasets}
