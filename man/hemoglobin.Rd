\name{hemoglobin}
\alias{hemoglobin}
\docType{data}
\title{
Treatment of insulin dependent diabetic children
}
\description{
16 insulin-dependent diabetic children were enrolled in a study involving
a new treatment. 8 children received the new treatment(N) while the
other 8 received the standard treatment(S).  The age and sex of the
child was recorded along with the measured value of gycosolated
hemoglobin both before and after treatment.
}
\usage{data(hemoglobin)}
\format{
  A data frame with 16 observations on the following 5 variables.
  \describe{
    \item{\code{age}}{age in years}
    \item{\code{sex}}{a factor with levels \code{F} \code{M}}
    \item{\code{treatment}}{a factor with levels \code{N} \code{S}}
    \item{\code{pre}}{measured value of hemoglobin before treatment}
    \item{\code{post}}{measured value of hemoglobin after treatment}
  }
}
\source{
Unknown
}
\examples{
data(hemoglobin)
## maybe str(hemoglobin) ; plot(hemoglobin) ...
}
\keyword{datasets}
