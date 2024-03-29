% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/faraway.R
\name{vif}
\alias{vif}
\alias{vif.default}
\alias{vif.lm}
\title{vif}
\usage{
vif(object)

\method{vif}{default}(object)

\method{vif}{lm}(object)
}
\arguments{
\item{object}{a data matrix (design matrix without intercept) or a model
object}
}
\value{
variance inflation factors
}
\description{
vif
}
\details{
Computes the variance inflation factors
}
\examples{

data(stackloss)
vif(stackloss[,-4])
#  Air.Flow Water.Temp Acid.Conc.
#    2.9065     2.5726     1.3336

}
\author{
Julian Faraway
}
