\name{exa}
\alias{exa}
\docType{data}
\title{
Simulated non-parametric regression data
}
\description{
True function is f(x)=sin^3(2pi x^3).
}
\usage{data(exa)}
\format{
  A data frame with 256 observations on the following 3 variables.
  \describe{
    \item{\code{x}}{input}
    \item{\code{y}}{response}
    \item{\code{m}}{true value}
  }
}
\source{
Haerdle, W. (1991). Smoothing Techniques with Implementation in S. New York:Springer.
}
\examples{
data(exa)
## maybe str(exa) ; plot(exa) ...
}
\keyword{datasets}
