\name{amlxray}
\alias{amlxray}
\docType{data}
\title{Match pair study for AML and Xray link}
\description{
  A matched case control study carried out to investigate the connection
  between X-ray usage and acute myeloid leukemia in childhood. The pairs
  are matched by age, race and county of residence.
}
\usage{data(amlxray)}
\format{
  A data frame with 238 observations on the following 11 variables.
  \describe{
    \item{\code{ID}}{a factor denoting the matched pairs}
    \item{\code{disease}}{0=control, 1=case}
    \item{\code{Sex}}{ \code{F} or \code{M}}
    \item{\code{downs}}{Presence of Downs syndrome: \code{no} or \code{yes}}
    \item{\code{age}}{Age in years}
    \item{\code{Mray}}{Did the mother ever have an Xray: \code{no} or \code{yes}}
    \item{\code{MupRay}}{Did the mother have an Xray of the upper body
      during pregnancy: \code{no} or \code{yes}}
    \item{\code{MlowRay}}{Did the mother have an Xray of the lower body
      during pregnancy: \code{no} or \code{yes}}
    \item{\code{Fray}}{Did the father ever have an Xray: \code{no} or \code{yes}}
    \item{\code{Cray}}{Did the child ever have an Xray: \code{no} or \code{yes}}
    \item{\code{CnRay}}{Total number of Xrays of the child \code{1}=none
      < \code{2}=1 or 2 < \code{3}=3 or 4 < \code{4}= 5 or more}
  }
}
\source{
  Chap T. Le (1998) "Applied Categorical Data Analysis" Wiley.
}
\keyword{datasets}
