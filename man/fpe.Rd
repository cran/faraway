\name{fpe}
\alias{fpe}
\title{1981 French Presidential Election}
\description{Elections for the French presidency proceed in two rounds. In 1981,
there were 10 candidates in the first round. The top two candidates
then went on to the second round, which was won by Francois
Mitterand over Valery Giscard-d'Estaing. The losers in the first
round can gain political favors by urging their supporters to vote for
one of the two finalists. Since voting is private, we cannot know how
these votes were transferred, we might hope to infer from the
published vote totals how this might have happened.  Data is given
for vote totals in every fourth department of France:
}
\usage{data(fpe)
}
\format{This dataframe contains the following columns (vote totals
  are in thousands)
  \describe{
    \item{\code{EI}}{ Electeur Inscrits (registered voters)}

    \item{\code{A}}{ Voters for Mitterand in the first round}
    \item{\code{B}}{ Voters for Giscard in the first round}
    \item{\code{C}}{ Voters for Chirac in the first round}
    \item{\code{D}}{ Voters for Communists in the first round}
    \item{\code{E}}{ Voters for Ecology party in the first round}
    \item{\code{F}}{ Voters for party F in the first round}
    \item{\code{G}}{ Voters for party G in the first round}
    \item{\code{H}}{ Voters for party H in the first round}
    \item{\code{I}}{ Voters for party I in the first round}
    \item{\code{J}}{ Voters for party J in the first round}
    \item{\code{K}}{ Voters for party K in the first round}
    \item{\code{A2}}{ Voters for Mitterand in the second round}
    \item{\code{B2}}{ Voters for party Giscard in the second round}
    \item{\code{N}}{ Difference between the number of voters in the second
      round and in the first round}

  }
}

\source{"The Teaching of Practical Statistics" by C.W. Anderson and
  R.M. Loynes, Wiley,1987}
\keyword{datasets}
