\name{logit}
\alias{logit}
\title{Logit transformation}
\description{Computes the logit transformation
}
\usage{logit(x)
}
\arguments{
  \item{x}{a numeric vector}
}
\value{log(x/(1-x))
}
\details{x <=0 or >=1 will return NA
}
\references{
}
\seealso{ilogit
}
\examples{logit(c(0.1,0.5,1.0,1.1))
#[1] -2.197225  0.000000        NA        NA
}
\author{Julian Faraway}
\keyword{manip}
