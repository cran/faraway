\name{suicide}
\alias{suicide}
\docType{data}
\title{Suicide method data from the UK}
\description{
  One year of suicide data from the United Kingdom crossclassified by
  sex, age and method.
}
\usage{data(suicide)}
\format{
  A data frame with 36 observations on the following 4 variables.
  \describe{
    \item{\code{y}}{number of people}
    \item{\code{cause}}{method used - a factor with levels \code{drug}
      (suicide by solid or liquid matter), \code{gas}, \code{gun} (guns,
      knives or explosives)
      \code{hang} (hanging, strangling, suffocating or drowning, \code{jump} \code{other}}
    \item{\code{age}}{a factor with levels \code{m} (middle-aged)
      \code{o} (old) \code{y} (young) }
    \item{\code{sex}}{a factor with levels \code{f} \code{m}}
  }
}
\source{
Everitt B. & Dunn G. (1991) "Applied Multivariate Data Analysis" Edward Arnold
}
\keyword{datasets}
