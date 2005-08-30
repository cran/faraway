\name{hprice}
\alias{hprice}
\docType{data}
\title{Housing prices in US cities 86-94}
\description{
  Data on housing prices in 36 US metropolitan statistical areas (MSAs)
over 9 years from 1986-1994 were collected.
}
\usage{data(hprice)}
\format{
  A data frame with 324 observations on the following 8 variables.
  \describe{
    \item{\code{narsp}}{log average sale price}
    \item{\code{ypc}}{average per capita income}
    \item{\code{perypc}}{percentage growth in per capita income}
    \item{\code{regtest}}{Regulatory environment index (high values = more regulations)}
    \item{\code{rcdum}}{Rent control - a factor with levels \code{0}=no \code{1}=yes}
    \item{\code{ajwtr}}{Adjacent to a coastline - a factor with levels \code{0}=no \code{1}=yes}
    \item{\code{msa}}{indicator for the MSA}
    \item{\code{time}}{Year 1=1986 to 9=1994}
  }
}
\source{
Longitudinal and Panel Data: Analysis and Applications in the Social Sciences, by Edward W. Frees, Cambridge University Press, August 2004.
}
\keyword{datasets}
