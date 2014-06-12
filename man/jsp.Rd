\name{jsp}
\alias{jsp}
\docType{data}
\title{
Junior School Project
}
\description{
Junior School Project collected
from primary (U.S. term is elementary) schools in inner London.
}
\usage{data(jsp)}
\format{
  A data frame with 3236 observations on the following 9 variables.
  \describe{
    \item{\code{school}}{50 schools code 1-50}
    \item{\code{class}}{a factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{gender}}{a factor with levels \code{boy} \code{girl}}
    \item{\code{social}}{class
of the father I=1; II=2; III nonmanual=3; III manual=4; IV=5; V=6;
Long-term unemployed=7; Not currently employed=8; Father absent=9}
    \item{\code{raven}}{test score}
    \item{\code{id}}{student id coded 1-1402}
    \item{\code{english}}{score on English}
    \item{\code{math}}{score on Maths}
    \item{\code{year}}{year of school}
  }
}
\source{
Mortimore, P., P. Sammons, L. Stoll, D. Lewis, and R. Ecob (1988). School Matters. Wells, UK: Open Books.
}
\references{
Goldstein, H. (1995). Multilevel Statistical Models (2 ed.). London: Arnold.
}
\examples{
data(jsp)
## maybe str(jsp) ; plot(jsp) ...
}
\keyword{datasets}
