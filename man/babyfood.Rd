\name{babyfood}
\alias{babyfood}
\docType{data}
\title{
Respiratory disease rates of babies fed in different ways
}
\description{
Study on
infant respiratory disease, namely the proportions of children
developing bronchitis or pneumonia in their first year of life by type
of feeding and sex.
}
\usage{data(babyfood)}
\format{
  A data frame with 6 observations on the following 4 variables.
  \describe{
    \item{\code{disease}}{number with disease}
    \item{\code{nondisease}}{number without disease}
    \item{\code{sex}}{a factor with levels \code{Boy} \code{Girl}}
    \item{\code{food}}{a factor with levels \code{Bottle} \code{Breast} \code{Suppl}}
  }
}
\source{
Payne, C. (1987). The GLIM System Release 3.77 Manual (2 ed.). Oxford: Nu- merical Algorithms Group.
}
\examples{
data(babyfood)
## maybe str(babyfood) ; plot(babyfood) ...
}
\keyword{datasets}
