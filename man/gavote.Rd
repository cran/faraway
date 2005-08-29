\name{gavote}
\alias{gavote}
\docType{data}
\title{Undercounted votes in Georgia in 2000 presidential election }
\description{
  The data comes from the US presidential election in the state of
  Georgia. The undercount is the difference between the number of ballots cast and votes
recorded. Voters may have chosen not to vote for president, voted
for more than one candidate (disqualified) or the equipment may have failed to register
their choice.
}
\usage{data(gavote)}
\format{
  A data frame with 159 observations on the following 10 variables. Each
  case represents a county in Georgia.
  \describe{
    \item{\code{equip}}{The voting equipment used: \code{LEVER},
      \code{OS-CC} (optical, central count), \code{OS-PC} (optical,
      precinct count) \code{PAPER}, \code{PUNCH}}
    \item{\code{econ}}{economic status of county: \code{middle} \code{poor} \code{rich}}
    \item{\code{perAA}}{percent of African Americans in county}
    \item{\code{rural}}{indicator of whether county is \code{rural} or \code{urban}}
    \item{\code{atlanta}}{indicator of whether county is in
      \code{Atlanta} or not: \code{notAtlanta}}
    \item{\code{gore}}{number of votes for Gore}
    \item{\code{bush}}{number of votes for Bush}
    \item{\code{other}}{number of votes for other candidates}
    \item{\code{votes}}{number of votes}
    \item{\code{ballots}}{number of ballots}
  }
}
\source{
   Meyer M. (2002) Uncounted Votes: Does Voting Equipment Matter?  Chance, 15(4), 33-38
}
\keyword{datasets}
