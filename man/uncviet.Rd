\name{uncviet}
\alias{uncviet}
\docType{data}
\title{UNC student opinions about the Vietnam War}
\description{
  A student newspaper conducted a survey of student opinions about the
  Vietnam War in May 1967. Responses were classified by sex, year in the
  program and one of four opinions. The survey was voluntary.
}
\usage{data(uncviet)}
\format{
  A data frame with 40 observations on the following 4 variables.
  \describe{
    \item{\code{y}}{the count}
    \item{\code{policy}}{a factor with levels
      \code{A} (defeat power of North Vietnam by widespread bombing
      and land invasion) \code{B} (follow the present policy)
      \code{C} (withdraw troops to strong points and open
   negotiations on elections involving the Viet Cong) \code{D} (immediate withdrawal of all U.S. troops)}
    \item{\code{sex}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{year}}{a factor with levels \code{Fresh} \code{Grad} \code{Junior} \code{Senior} \code{Soph}}
  }
}
\source{
  M. Aitkin and D. Anderson and B. Francis and J. Hinde
  (1989) "Statistical Modelling in GLIM" Oxford University Press.
}
\keyword{datasets}
