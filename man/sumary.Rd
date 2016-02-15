\name{sumary}
\alias{sumary}
\alias{sumary,lm-method}
\alias{sumary,glm-method}
\alias{sumary,merMod-method}
\title{Abbreviated Regression Summary}
\description{Generic summaries for lm, glm and mer objects}
\usage{
sumary(object, ...)
}
\arguments{
  \item{object}{An lm, glm or mer object returned from lm(), glm() or
    lmer() respectively}
  \item{...}{further arguments passed to or from other methods.}
}


\value{returns the same as \code{summary()}
}

\details{This generic function provides an abbreviated regression output
  containing the more useful information. Users wanting to see more
  are advised to use \code{summary()}
}


\references{This function is adapted from the \code{display()} function in the
  \code{arm} package}

\seealso{\code{\link[base]{summary}},
         \code{\link[stats]{lm}},
         \code{\link[stats]{glm}},
         \code{\link[lme4]{lmer}}
       }

\examples{
data(stackloss)
object <- lm(stack.loss ~ .,stackloss)
sumary(object)
}
\author{Julian Faraway}
\keyword{regression}


