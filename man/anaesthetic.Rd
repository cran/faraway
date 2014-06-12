\name{anaesthetic}
\alias{anaesthetic}
\docType{data}
\title{
Time in minutes to eye opening after reversal of anaesthetic.
}
\description{
A doctor at major London hospital compared the
effects of 4 anaesthetics used in major operations.
80 patients were divided into groups of 20.
}
\usage{data(anaesthetic)}
\format{
  A data frame with 80 observations on the following 2 variables.
  \describe{
    \item{\code{breath}}{time in minutes to start breathing unassisted}
    \item{\code{tgrp}}{Four treatment groups \code{A} \code{B} \code{C} \code{D}}
  }
}
\source{
Chatfield C. (1995) Problem Solving: A Statistician's Guide, 2ed Chapman Hall.
}
\examples{
data(anaesthetic)
## maybe str(anaesthetic) ; plot(anaesthetic) ...
}
\keyword{datasets}
