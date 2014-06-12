\name{crawl}
\alias{crawl}
\docType{data}
\title{
Crawling babies by month
}
\description{
A study investigated whether babies take longer to learn to
  crawl in cold months when they are often bundled in clothes that
  restrict their movement, than in warmer months. The study sought an
  association between babies' first crawling age and the average
  temperature during the month they first try to crawl (about 6 months
  after birth). Parents brought their babies into the University of
  Denver Infant Study Center between 1988-1991 for the study. The
  parents reported the birth month and age at which their child was
  first able to creep or crawl a distance of four feet in one minute.
  Data were collected on 208 boys and 206 girls (40 pairs of which
  were twins)
}
\usage{data(crawl)}
\format{
  A data frame with 12 observations on the following 4 variables.
  \describe{
    \item{\code{crawling}}{average crawling age in weeks}
    \item{\code{SD}}{standard deviation of crawling age}
    \item{\code{n}}{sample size}
    \item{\code{temperature}}{average temperature(F) six months after birth}
  }
}
\source{
Benson, Janette. (1993). Infant Behavior and Development
}
\examples{
data(crawl)
## maybe str(crawl) ; plot(crawl) ...
}
\keyword{datasets}
