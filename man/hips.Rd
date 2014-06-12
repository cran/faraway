\name{hips}
\alias{hips}
\docType{data}
\title{
Ankylosing Spondylitis
}
\description{
Data from Royal Mineral Hospital in Bath. AS is a chronic form of arthritis. A study conducted  to
  determine whether daily stretching of the hip tissues would improve
  mobility.  39 ``typical'' AS patients were randomly allocated to
  control (standard treatment) group or the treatment group in a 1:2
  ratio.  Responses were flexion and rotation angles at the hip
  measured in degrees. Larger numbers indicate more flexibility.
}
\usage{data(hips)}
\format{
  A data frame with 78 observations on the following 7 variables.
  \describe{
    \item{\code{fbef}}{flexion angle before}
    \item{\code{faft}}{flexion angle after}
    \item{\code{rbef}}{rotation angle before}
    \item{\code{raft}}{rotation angle after}
    \item{\code{grp}}{treatment group - a factor with levels \code{control} \code{treat}}
    \item{\code{side}}{side of the body - a factor with levels \code{right} \code{left}}
    \item{\code{person}}{id for the individual}
  }
}
\source{
Chatfield C. (1995) Problem Solving: A Statistician's Guide, 2ed Chapman Hall.
}
\examples{
data(hips)
## maybe str(hips) ; plot(hips) ...
}
\keyword{datasets}
