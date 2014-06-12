\name{gammaray}
\alias{gammaray}
\docType{data}
\title{
Xray decay from a gamma ray burst
}
\description{
The X-ray decay light curve of Gamma ray burst 050525a obtained with the X-Ray Telescope (XRT) on board the Swift satellite. The dataset has 63 brightness measurements in the 0.4-4.5 keV spectral band at times ranging from 2 minutes to 5 days after the burst.
}
\usage{data(gammaray)}
\format{
  A data frame with 63 observations on the following 3 variables.
  \describe{
    \item{\code{time}}{in seconds since burst}
    \item{\code{flux}}{X-ray flux in units of 10^{-11} erg/cm2/s, 2-10 keV}
    \item{\code{error}}{measurement error of the flux based on detector signal-to-noise values}
  }
}
\source{
A. J. Blustin and 64 coauthors, Astrophys. J. 637, 901-913 2006. Available at http://arxiv.org/abs/astro-ph/0507515.
}
\examples{
data(gammaray)
## maybe str(gammaray) ; plot(gammaray) ...
}
\keyword{datasets}
