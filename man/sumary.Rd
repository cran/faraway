\name{sumary}
\alias{sumary}
\title{Abbreviated Regression Summary}
\description{Abbreviated Regression Summary intended to replace base summary
}
\usage{
sumary(object)
}
\arguments{

  \item{object}{An lm object returned from lm()}

}
\value{none
}
\examples{
data(stackloss)
object <- lm(stack.loss ~ .,stackloss)
sumary(object)
}
\author{Julian Faraway}
\keyword{regression}
