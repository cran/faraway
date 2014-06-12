\name{airpass}
\alias{airpass}
\docType{data}
\title{
Airline passengers
}
\description{
Monthly totals of airline passengers from 1949 to 1951
}
\usage{data(airpass)}
\format{
  A data frame with 144 observations on the following 2 variables.
  \describe{
    \item{\code{pass}}{number of passengers in thousands}
    \item{\code{year}}{the date as a decimal}
  }
}
\details{
Well known time series example dataset
}
\source{
Brown, R.G.(1962)
Smoothing, Forecasting and Prediction of Discrete Time Series.
 Englewood Cliffs, N.J.: Prentice-Hall.
}
\references{
Box, G.E.P., Jenkins, G.M. and Reinsel, G.C. (1994)
Time Series Analysis, Forecasting and Control, 3rd edn.
 Englewood Cliffs, N.J.: Prentice-Hall.
}
\examples{
data(airpass)
## maybe str(airpass) ; plot(airpass) ...
}
\keyword{datasets}
