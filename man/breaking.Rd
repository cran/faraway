\name{breaking}
\alias{breaking}
\docType{data}
\title{Breaking strength of materials}
\description{
  An experiment was conducted to select the supplier of raw materials
  for production of a component. The breaking strength of the component
  was the objective of interest.  Four suppliers were considered. The
  four operators can only produce one component each per day. A Latin
  square design was used.
}
\usage{data(breaking)}
\format{
  A data frame with 16 observations on the following 4 variables.
  \describe{
    \item{\code{y}}{The breaking strength of the component}
    \item{\code{operator}}{the operator - a factor with levels \code{op1} \code{op2} \code{op3} \code{op4}}
    \item{\code{day}}{the day of production - a factor with levels \code{day1} \code{day2} \code{day3} \code{day4}}
    \item{\code{supplier}}{the supplier of the raw material - a factor with levels \code{A} \code{B} \code{C} \code{D}}
  }
}
\source{
  Lentner M. and Bishop T. (1986) Experimental Design and Analysis,
  Valley Book Company
}
\keyword{datasets}
