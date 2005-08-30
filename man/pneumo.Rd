\name{pneumo}
\alias{pneumo}
\docType{data}
\title{Pneumonoconiosis in coal miners}
\description{
The data for this example contains the number of coal miners
classified by radiological examination into one of three categories of
pneumonoultramicroscopicosilicovolcanoconiosis (known as
pneumonoconiosis for short) and by number of years spent
working at the coal face divided into eight categories. }
\usage{data(pneumo)}
\format{
  A data frame with 24 observations on the following 3 variables.
  \describe{
    \item{\code{Freq}}{number of miners}
    \item{\code{status}}{pneumoconiosis status - a factor with levels \code{mild} \code{normal} \code{severe}}
    \item{\code{year}}{number of years service (midpoint of interval)}
  }
}
\source{
  M. Aitkin and D. Anderson and B. Francis and J. Hinde
  (1989) "Statistical Modelling in GLIM" Oxford University Press.
}
\keyword{datasets}
