\name{abrasion}
\alias{abrasion}
\title{Wear on materials according to type, run and position}
\description{The \code{composite} data frame has 16 rows and 4 columns.
  Four materials were fed into a wear testing machine
  and the amount of wear recorded. Four samples could be processed
  at the same time and the position of these samples may be important.
  A Latin square design was used.
}
\usage{
data(abrasion)
}
\format{
  This data frame contains the following columns:
  \describe{

    \item{\code{run}}{
      The run number 1-4
    }

    \item{\code{position}}{
      The position number 1-4
    }

    \item{\code{material}}{
      The material A-D
    }

    \item{\code{wear}}{
      The wear measured loss of weight in 0.1mm over testing period
    }}}
\source{The Design and Analysis of Industrial Experiments by O. Davies,
  1954, published by Wiley}
  
\keyword{datasets}

