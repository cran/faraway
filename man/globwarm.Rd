\name{globwarm}
\alias{globwarm}
\docType{data}
\title{
Northern Hemisphere temperatures and climate proxies in the last millenia
}
\description{
Average Northen Hemisphere Temperature from 1856-2000 and eight climate
proxies from 1000-2000AD. Data can be used to predict temperatures prior
to 1856.
}
\usage{data(globwarm)}
\format{
  A data frame with 1001 observations on the following 10 variables.
  \describe{
    \item{\code{nhtemp}}{Northern hemisphere average temperature (C) provided by the UK Met Office (known as HadCRUT2)}
    \item{\code{wusa}}{Tree ring proxy information from the Western USA.}
    \item{\code{jasper}}{Tree ring proxy information from Canada.}
    \item{\code{westgreen}}{Ice core proxy information from west Greenland}
    \item{\code{chesapeake}}{Sea shell proxy information from Chesapeake Bay}
    \item{\code{tornetrask}}{Tree ring proxy information from Sweden}
    \item{\code{urals}}{Tree ring proxy information from the Urals}
    \item{\code{mongolia}}{Tree ring proxy information from Mongolia}
    \item{\code{tasman}}{Tree ring proxy information from Tasmania}
    \item{\code{year}}{Year 1000-2000AD}
  }
}
\details{
See the source and references below for the original data. Only some
proxies have been included here. Some missing values have been imputed.
The proxy data have been smoothed. This version of the data is intended
only for demonstration purposes. If you are specifically interested in
the subject matter, use the original data.
}
\source{
P.D. Jones and M.E. Mann (2004) "Climate Over Past Millennia" Reviews of Geophysics, Vol. 42, No. 2, RG2002, doi:10.1029/2003RG000143
}
\references{
www.ncdc.noaa.gov/paleo/pubs/jones2004/jones2004.html
}
\examples{
data(globwarm)
## maybe str(globwarm) ; plot(globwarm) ...
}
\keyword{datasets}
