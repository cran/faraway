\name{vif}
\alias{vif}
\alias{vif.default}
\alias{vif.lm}
\title{vif}
\description{Computes the variance inflation factors
}
\usage{vif(object)
}
\arguments{
  \item{object}{a data matrix (design matrix without intercept) or a model object}
}
\value{variance inflation factors
}
\examples{
data(stackloss)
vif(stackloss[,-4])
#  Air.Flow Water.Temp Acid.Conc. 
#    2.9065     2.5726     1.3336 
}
\author{Julian Faraway}
