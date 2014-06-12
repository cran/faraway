\name{newhamp}
\alias{newhamp}
\docType{data}
\title{
New Hampshire Democratic Party Primary 2008
}
\description{
Votes and other demographic information from 276 wards in the 2008
Democratic Party presidential primary.
}
\usage{data(newhamp)}
\format{
  A data frame with 276 observations on the following 12 variables.
  \describe{
    \item{\code{votesys}}{The voting system used where H is counted by hand
      and D is counted by machine.}
    \item{\code{Obama}}{The number of votes for Barack Obama.}
    \item{\code{Clinton}}{The number of votes for Hillary Clinton.}
    \item{\code{dem}}{The total number of votes cast in the Democratic primary (there were
  other candidates besides Clinton and Obama).}
      \item{\code{povrate}}{The poverty rate as a proportion as determined by the 2000 census.}
    \item{\code{pci}}{Per capita annual income in USD in 1999.}
    \item{\code{Dean}}{The proportion of voters for Howard Dean in the 2004 Democratic primary.}
    \item{\code{Kerry}}{The proportion of voters for John Kerry in the 2004 Democratic primary.}
    \item{\code{white}}{The proportion of non-Hispanic whites according to the 2000 census.}
    \item{\code{absentee}}{The proportion voting by absentee ballot.}
    \item{\code{population}}{An estimate of the population from 2002.}
    \item{\code{pObama}}{Proportion voting for Obama}
  }
}
\details{
On the 8th January 2008, primaries to select US presidential candidates
were held in New Hampshire. In the Democratic party primary, Hillary Clinton
defeated Barack Obama contrary to the expectations pre-election opinion
polls. Essentially two different voting technologies were used in New Hampshire.
Some wards used paper ballots, counted by hand while others used optically
scanned ballots, counted by machine. Among the paper ballots, Obama had more
votes than Clinton while Clinton defeated Obama on just the machine counted
ballots. Since the method of voting should make no causal difference to
the outcome, suspicions have been raised regarding the integrity of the
election. 
}
\source{
Herron, M., W. M. Jr, and J. Wand (2008). Voting Technology and the 2008 New Hampshire Primary. Wm. & Mary Bill Rts. J. 17, 351-374.
}
\keyword{datasets}
