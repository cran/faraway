\name{fortune}
\alias{fortune}
\docType{data}
\title{
Billionaires' wealth and age
}
\description{
Fortune magazine publishes a list of the world's billionaires
  each year. The 1992 list includes 233 individuals. Their wealth,
  age, and geographic location (Asia, Europe, Middle East, United
  States, and Other) are reported. 
}
\usage{data(fortune)}
\format{
  A data frame with 232 observations on the following 3 variables.
  \describe{
    \item{\code{wealth}}{Billions of dollars}
    \item{\code{age}}{age in years}
    \item{\code{region}}{a factor with levels \code{A}, Asia, \code{E}, Europe, \code{M}, Middle East, \code{O} Other, \code{U} USA}
  }
}
\source{
Fortune magazine
}
\examples{
data(fortune)
## maybe str(fortune) ; plot(fortune) ...
}
\keyword{datasets}
