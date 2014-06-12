\name{vision}
\alias{vision}
\docType{data}
\title{
Acuity of vision in response to light flash
}
\description{
The acuity of vision for seven subjects was tested. The response
is the lag in milliseconds between a light flash and a response in
the cortex of the eye. Each eye is tested at four different powers of lens.
An object at the distance of the second number appears to be at distance
of the first number.
}
\usage{data(vision)}
\format{
  A data frame with 56 observations on the following 4 variables.
  \describe{
    \item{\code{acuity}}{a numeric vector}
    \item{\code{power}}{a factor with levels \code{6/6} \code{6/18} \code{6/36} \code{6/60}}
    \item{\code{eye}}{a factor with levels \code{left} \code{right}}
    \item{\code{subject}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7}}
  }
}
\source{
Crowder, M. J. and D. J. Hand (1990). Analysis of Repeated Measures. London:
Chapman & Hall.
}
\examples{
data(vision)
## maybe str(vision) ; plot(vision) ...
}
\keyword{datasets}
