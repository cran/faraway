\name{troutegg}
\alias{troutegg}
\docType{data}
\title{Survival of trout eggs depending on time and location}
\description{
  Boxes of trout eggs were buried at five different stream locations
  and retrieved at 4 different times. The number of surviving eggs was
  recorded. The box was not returned to the stream.
}
\usage{data(troutegg)}
\format{
  A data frame with 20 observations on the following 4 variables.
  \describe{
    \item{\code{survive}}{the number of surviving eggs}
    \item{\code{total}}{the number of eggs in the box}
    \item{\code{location}}{the location in the stream with levels \code{1} \code{2} \code{3} \code{4} \code{5}}
    \item{\code{period}}{the number of weeks after placement that the
      box was withdrawn levels \code{4} \code{7} \code{8} \code{11}}
  }
}
\source{
  Manly B. (1978) Regression models for proportions with extraneous
  variance. Biometrie-Praximetrie, 18, 1-18.
}
\references{
  Hinde J. and Demetrio C. (1988) Overdispersion: Models and
  estimation. Computational Statistics and Data Analysis. 27, 151-170.
}
\keyword{datasets}
