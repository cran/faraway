\name{qqnorml}
\alias{qqnorml}
\title{Labeled QQ plot}
\description{Makes a labeled QQ plot
}
\usage{
qqnorml(y,main = "Normal Q-Q Plot", xlab = "Theoretical Quantiles", 
    ylab = "Sample Quantiles",...)
}
\arguments{

  \item{y}{A numeric vector}

  \item{main}{main label}

  \item{xlab}{x-axis label}
  
  \item{ylab}{y-axis label}

  \item{...}{arguments passed to plot()}
}
\value{none
}
\details{
}
\references{
}
\seealso{
}
\examples{
qqnorml(rnorm(16))
}
\author{Julian Faraway}
\keyword{}
