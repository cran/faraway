\name{ilogit}
\alias{ilogit}
\title{Inverse Logit Transformation}
\description{
  Computes the inverse logit transformation
}
\usage{ilogit(x)
}
\arguments{
  \item{x}{a numeric vector
}}
\value{
  exp(x)/(1+exp(x))
}
\seealso{logit
}
\examples{
ilogit(1:3)
#[1] 0.7310586 0.8807971 0.9525741
}
\author{Julian Faraway}

