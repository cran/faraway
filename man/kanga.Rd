\name{kanga}
\alias{kanga}
\docType{data}
\title{
Kangaroo skull measurements
}
\description{
Sex and species of an specimens of
kangaroo.
}
\usage{data(kanga)}
\format{
  A data frame with 148 observations on the following 20 variables.
  \describe{
    \item{\code{species}}{a factor with levels \code{fuliginosus} \code{giganteus} \code{melanops}}
    \item{\code{sex}}{a factor with levels \code{Female} \code{Male}}
    \item{\code{basilar.length}}{a numeric vector}
    \item{\code{occipitonasal.length}}{a numeric vector}
    \item{\code{palate.length}}{a numeric vector}
    \item{\code{palate.width}}{a numeric vector}
    \item{\code{nasal.length}}{a numeric vector}
    \item{\code{nasal.width}}{a numeric vector}
    \item{\code{squamosal.depth}}{a numeric vector}
    \item{\code{lacrymal.width}}{a numeric vector}
    \item{\code{zygomatic.width}}{a numeric vector}
    \item{\code{orbital.width}}{a numeric vector}
    \item{\code{.rostral.width}}{a numeric vector}
    \item{\code{occipital.depth}}{a numeric vector}
    \item{\code{crest.width}}{a numeric vector}
    \item{\code{foramina.length}}{a numeric vector}
    \item{\code{mandible.length}}{a numeric vector}
    \item{\code{mandible.width}}{a numeric vector}
    \item{\code{mandible.depth}}{a numeric vector}
    \item{\code{ramus.height}}{a numeric vector}
  }
}
\source{
Andrews and Herzberg (1985) Chapter 53.
}
\references{
Andrews, D. F. and Herzberg, A. M. (1985). Data. Springer-Verlag, New York.
}
\examples{
data(kanga)
## maybe str(kanga) ; plot(kanga) ...
}
\keyword{datasets}
