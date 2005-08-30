\name{broccoli}
\alias{broccoli}
\docType{data}
\title{Broccoli weight variation}
\description{
  A number of growers supply broccoli to a food processing plant. The
plant instructs the growers to pack the broccoli into standard size boxes.
There should be 18 clusters of broccoli per box and each cluster should
weigh between 1.33 and 1.5 pounds. Because the growers use different
varieties, methods of cultivation etc, there is some variation in the
cluster weights. The plant manager selected 3 growers at random and then
4 boxes at random supplied by these growers. 3 clusters were selected
from each box.
}
\usage{data(broccoli)}
\format{
  A data frame with 36 observations on the following 4 variables.
  \describe{
    \item{\code{wt}}{weight of broccoli}
    \item{\code{grower}}{the grower - a factor with levels \code{1} \code{2} \code{3}}
    \item{\code{box}}{the box - a factor with levels \code{1} \code{2} \code{3} \code{4}}
    \item{\code{cluster}}{the cluster - a factor with levels \code{1} \code{2} \code{3}}
  }
}
\source{
  Lentner M. and Bishop T. (1986) Experimental Design and Analysis,
  Valley Book Company
}
\keyword{datasets}
