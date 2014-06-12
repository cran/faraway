\name{worldcup}
\alias{worldcup}
\docType{data}
\title{
Data on players from the 2010 World Cup
}
\description{
Data on players from the 2010 World Cup
}
\usage{data(worldcup)}
\format{
  A data frame with 595 observations on the following 7 variables.
  \describe{
    \item{\code{Team}}{Country}
    \item{\code{Position}}{a factor with levels \code{Defender} \code{Forward} \code{Goalkeeper} \code{Midfielder}}
    \item{\code{Time}}{Time played in minutes}
    \item{\code{Shots}}{Number of shots attempted}
    \item{\code{Passes}}{Number of passes made}
    \item{\code{Tackles}}{Number of tackles made}
    \item{\code{Saves}}{Number of saves made}
  }
}
\details{
None
}
\source{
Lost
}
\examples{
data(worldcup)
## maybe str(worldcup) ; plot(worldcup) ...
}
\keyword{datasets}
