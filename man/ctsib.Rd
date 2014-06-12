\name{ctsib}
\alias{ctsib}
\docType{data}
\title{
Effects of surface and vision on balance.
}
\description{
An experiment was conducted to study the effects of surface and vision
on balance. The balance of subjects were observed for two different
surfaces and for restricted and unrestricted vision. Balance was
assessed qualitatively on an ordinal four-point scale based on
observation by the experimenter. Forty subjects were studied, twenty males
and twenty females ranging in age from 18 to 38, with heights given
in cm and weights in kg. The subjects were tested while standing on
foam or a normal surface and with their eyes closed or open or with 
a dome placed over their head.  Each subject was tested twice in each
of the surface and eye combinations for a total of 12 measures per
subject.
}
\usage{data(ctsib)}
\format{
  A data frame with 480 observations on the following 8 variables.
  \describe{
    \item{\code{Subject}}{an indicator}
    \item{\code{Sex}}{a factor with levels \code{female} \code{male}}
    \item{\code{Age}}{in years}
    \item{\code{Height}}{in cm}
    \item{\code{Weight}}{in kg}
    \item{\code{Surface}}{a factor with levels \code{foam} \code{norm}}
    \item{\code{Vision}}{a factor with levels \code{closed} \code{dome} \code{open}}
    \item{\code{CTSIB}}{a four point scale measuring balance}
  }
}
\source{
Steele, R. (1998). Effect of surface and vision on balance. Ph. D. thesis, Depart- ment of Physiotherapy, University of Queensland.
}
\references{
OzDasl
}
\examples{
data(ctsib)
## maybe str(ctsib) ; plot(ctsib) ...
}
\keyword{datasets}
