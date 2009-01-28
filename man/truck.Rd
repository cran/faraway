\name{truck}
\alias{truck}
\docType{data}
\title{Truck leaf spring experiment}
\description{
  Data on an experiment
concerning the production of leaf springs for trucks.
  A \eqn{2^{5-1}} fractional factorial experiment with 3 replicates was
  carried out with objective of recommending production settings to
  achieve a free height as close as possible to 8 inches.
}
\usage{data(truck)}
\format{
  A data frame with 48 observations on the following 6 variables.
  \describe{
    \item{\code{B}}{furnace temperature - a factor with levels \code{+} \code{-}}
    \item{\code{C}}{heating time - a factor with levels \code{+} \code{-}}
    \item{\code{D}}{transfer time - a factor with levels \code{+} \code{-}}
    \item{\code{E}}{hold-down time - a factor with levels \code{+} \code{-}}
    \item{\code{O}}{quench oil temperature - a factor with levels \code{+} \code{-}}
    \item{\code{height}}{leaf spring free height in inches}
  }
}
\source{
J. J. Pignatiello and J. S. Ramberg (1985) Contribution to discussion of
offline quality control, parameter design and the Taguchi method,
Journal of Quality Technology, \bold{17} 198-206.
}
\references{
P. McCullagh and J. Nelder (1989) "Generalized Linear Models" Chapman
and Hall, 2nd ed.
}
\keyword{datasets}
