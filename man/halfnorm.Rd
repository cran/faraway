\name{halfnorm}
\alias{halfnorm}
\title{Half Normal Plot}
\description{
  Makes a half-normal plot
}
\usage{halfnorm(x, nlab = 2, labs = as.character(1:length(x)), ylab = "Sorted Data", 
            ...) 
}
\arguments{
  \item{x}{a numeric vector
  }

  \item{nlab}{number of points to label
  }

  \item{labs}{labels for points
  }

  \item{ylab}{label for Y-axis
  }

  \item{...}{arguments passed to plot()}
}

\value{
  none
}
\details{
}
\references{
}
\seealso{
}
\examples{
halfnorm(runif(10))
}
\author{Julian Faraway}
\keyword{}
