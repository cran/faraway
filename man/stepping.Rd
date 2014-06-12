\name{stepping}
\alias{stepping}
\docType{data}
\title{
Stepping and effect on heart rate
}
\description{
An experiment was conducted
to explore the nature of the relationship between a person's heart
rate and the frequency at which that person stepped up and down on
steps of various heights.
}
\usage{data(stepping)}
\format{
  A data frame with 30 observations on the following 6 variables.
  \describe{
    \item{\code{Order}}{running order within the experiment}
    \item{\code{Block}}{Experimenter used}
    \item{\code{Height}}{0 if step at the low (5.75in) height, 1 if at the high (11.5in) height}
    \item{\code{Frequency}}{the rate of stepping. 0 if slow (14 steps/min), 1 if medium (21 steps/min), 2 if high (28 steps/min)}
    \item{\code{RestHR}}{the resting heart rate of the subject before a trial, in beats per minute}
    \item{\code{HR}}{ the final heart rate of the subject after a trial, in beats per minute }
  }
}
\source{
Unknown
}
\examples{
data(stepping)
## maybe str(stepping) ; plot(stepping) ...
}
\keyword{datasets}
