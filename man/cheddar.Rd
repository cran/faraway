\name{cheddar}
\alias{cheddar}
\docType{data}
\title{
Taste of Cheddar cheese
}
\description{
In a study of cheddar cheese from the
LaTrobe Valley of Victoria, Australia, samples of cheese were analyzed
for their chemical composition and were subjected to taste tests.
Overall taste scores were obtained by combining the scores from
several tasters.
}
\usage{data(cheddar)}
\format{
  A data frame with 30 observations on the following 4 variables.
  \describe{
    \item{\code{taste}}{a subjective taste score}
    \item{\code{Acetic}}{concentration of acetic acid (log scale)}
    \item{\code{H2S}}{concentration of hydrogen sulfice (log scale)}
    \item{\code{Lactic}}{concentration of lactic acid}
  }
}
\source{
Unknown
}
\examples{
data(cheddar)
## maybe str(cheddar) ; plot(cheddar) ...
}
\keyword{datasets}
