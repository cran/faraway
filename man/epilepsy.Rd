\name{epilepsy}
\alias{epilepsy}
\docType{data}
\title{
Epileptic seizures in clinical trial of drug
}
\description{
Data from a clinical trial
of 59 epileptics.  For a baseline, patients were observed for 8 weeks
and the number of seizures recorded. The patients were then randomized
to treatment by the drug Progabide (31 patients) or to the placebo
group (28 patients).  They were observed for four 2-week periods and
the number of seizures recorded.
}
\usage{data(epilepsy)}
\format{
  A data frame with 295 observations on the following 6 variables.
  \describe{
    \item{\code{seizures}}{number of seizures}
    \item{\code{id}}{identifying number}
    \item{\code{treat}}{1=treated, 0=not}
    \item{\code{expind}}{0=baseline period, 1=treatment period}
    \item{\code{timeadj}}{weeks of period}
    \item{\code{age}}{in years}
  }
}
\source{
Thall, P. F. and S. C. Vail (1990). Some covariance models for longitudinal count data with overdispersion. Biometrics 46, 657-671.
}
\references{
Breslow, N. E. and D. G. Clayton (1993). Approximate inference in generalized
linear mixed models. Journal of the American Statistical Association 88,
9-25.
Diggle, P. J., P. Heagerty, K. Y. Liang, and S. L. Zeger (2002). Analysis of Longitudinal Data (2 ed.). Oxford: Oxford University Press.
}
\examples{
data(epilepsy)
## maybe str(epilepsy) ; plot(epilepsy) ...
}
\keyword{datasets}
