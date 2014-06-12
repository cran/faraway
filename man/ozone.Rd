\name{ozone}
\alias{ozone}
\docType{data}
\title{
Ozone in LA in 1976
}
\description{
A study the relationship
between atmospheric ozone concentration and meteorology in the Los
Angeles Basin in 1976.  A number of cases with missing variables have
been removed for simplicity.
}
\usage{data(ozone)}
\format{
  A data frame with 330 observations on the following 10 variables.
  \describe{
    \item{\code{O3}}{Ozone conc., ppm, at Sandbug AFB.}
    \item{\code{vh}}{a numeric vector}
    \item{\code{wind}}{wind speed}
    \item{\code{humidity}}{a numeric vector}
    \item{\code{temp}}{temperature}
    \item{\code{ibh}}{inversion base height}
    \item{\code{dpg}}{Daggett pressure gradient}
    \item{\code{ibt}}{a numeric vector}
    \item{\code{vis}}{visibility}
    \item{\code{doy}}{day of the year}
  }
}
\source{
Breiman, L. and J. H. Friedman (1985). Estimating optimal transformations for
multiple regression and correlation. Journal of the American Statistical Association 80, 580-598.
}
\examples{
data(ozone)
## maybe str(ozone) ; plot(ozone) ...
}
\keyword{datasets}
