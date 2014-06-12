\name{phbirths}
\alias{phbirths}
\docType{data}
\title{
Birth weights in Philadelphia
}
\description{
Data based on a 5% sample of all births occurring in Philadelphia in 1990.
}
\usage{data(phbirths)}
\format{
  A data frame with 1115 observations on the following 5 variables.
  \describe{
    \item{\code{black}}{is the mother Black?}
    \item{\code{educ}}{mother's years of education}
    \item{\code{smoke}}{does the mother smoke during pregnancy?}
    \item{\code{gestate}}{gestational age in weeks}
    \item{\code{grams}}{birth weight in grams}
  }
}
\source{
I. T. Elo, G. Rodriguez and H. Lee (2001). Racial and Neighborhood Disparities in Birthweight in Philadelphia. Paper presented at the Annual Meeting of the Population Association of America, Washington, DC 2001.
}
\examples{
data(phbirths)
## maybe str(phbirths) ; plot(phbirths) ...
}
\keyword{datasets}
