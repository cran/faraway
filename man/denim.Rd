\name{denim}
\alias{denim}
\docType{data}
\title{
Denim wastage by supplier
}
\description{
Five suppliers cut denim material for a jeans manufacturer. An
  algorithm is used to estimate how much material will be wasted given
  the dimensions of the material supplied. Typically, a supplier wastes
  more material than the target based on the algorithm although occasionally
  they waste less. The percentage of waste relative to target was collected
  weekly for the 5 suppliers. In all, 95 observations were recorded. 
}
\usage{data(denim)}
\format{
  A data frame with 95 observations on the following 2 variables.
  \describe{
    \item{\code{waste}}{percentage wastage}
    \item{\code{supplier}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5}}
  }
}
\source{
Unknown
}
\examples{
data(denim)
## maybe str(denim) ; plot(denim) ...
}
\keyword{datasets}
