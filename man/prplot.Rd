\name{prplot}
\alias{prplot}
\title{Partial Residual Plot}
\description{Makes a Partial Residual plot
}
\usage{
prplot(g,i)
}
\arguments{

  \item{g}{An object returned from lm()}

  \item{i}{index of predictor}

}
\value{none
}
\examples{
data(stackloss)
g <- lm(stack.loss ~ .,stackloss)
prplot(g,1)
}
\author{Julian Faraway}
\keyword{regression}
