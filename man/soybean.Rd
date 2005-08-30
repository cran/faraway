\name{soybean}
\alias{soybean}
\docType{data}
\title{Germination failures for soybean seeds}
\description{
  An experiment was conducted to compare the germination rates of the
  five varieties of soybean. Five plots were available.
}
\usage{data(soybean)}
\format{
  A data frame with 25 observations on the following 3 variables.
  \describe{
    \item{\code{variety}}{the variety - a factor with levels \code{arasan} \code{check} \code{fermate} \code{semesan} \code{spergon}}
    \item{\code{replicate}}{the plot - a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5}}
    \item{\code{failure}}{the number of failures out of 100 planted seeds}
  }
}
\source{
  Snedecor G. and Cochran W. (1967) Statistical Methods (6th Ed) Iowa
  State University Press
}
\keyword{datasets}
