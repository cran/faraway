\name{wbca}
\alias{wbca}
\docType{data}
\title{Wisconsin breast cancer database}
\description{
  Data come from a study of breast cancer in Wisconsin. There
are 681 cases of potentially cancerous tumors of which 238 are actually
malignant. Determining whether a tumor is really malignant is traditionally
determined by an invasive surgical procedure. The purpose of this study was
to determine whether a new procedure called fine needle aspiration
which draws only a small sample of tissue could be effective in determining
tumor status.
}
\usage{data(wbca)}
\format{
  A data frame with 681 observations on the following 10 variables.
  \describe{
    \item{\code{Class}}{0 if malignant, 1 if benign}
    \item{\code{Adhes}}{marginal adhesion}
    \item{\code{BNucl}}{bare nuclei}
    \item{\code{Chrom}}{bland chromatin}
    \item{\code{Epith}}{epithelial cell size}
    \item{\code{Mitos}}{mitoses}
    \item{\code{NNucl}}{normal nucleoli}
    \item{\code{Thick}}{clump thickness}
    \item{\code{UShap}}{cell shape uniformity}
    \item{\code{USize}}{cell size uniformity}
  }
}
\details{
The predictor values are determined by a doctor observing the cells and
rating them on a scale from 1 (normal) to 10 (most abnormal) with respect
to the particular characteristic.
}
\source{
  Bennett, K.,P., and Mangasarian, O.L., Neural network training via linear programming. In P. M. Pardalos, editor, Advances in Optimization and Parallel Computing, pages 56-57. Elsevier Science, 1992
}
\keyword{datasets}
