\name{semicond}
\alias{semicond}
\non_function{}
\title{Semiconductor split-plot experiment}
\description{
    The \code{semicond} data frame has 48 rows and 5 columns.
}
\format{
    This data frame contains the following columns:
    \describe{
        \item{resistance}{
            a numeric vector
        }
        \item{ET}{
            a factor with levels \code{1} to \code{4} representing etch time.
        }
        \item{Wafer}{
            a factor with levels \code{1} to \code{3}
        }
        \item{position}{
            a factor with levels \code{1} to \code{4}
        }
        \item{Grp}{
            an ordered factor with levels
            \code{1/1} < \code{1/2} < \code{1/3} < \code{2/1} < \code{2/2} < \code{2/3} < \code{3/1} < \code{3/2} < \code{3/3} < \code{4/1} < \code{4/2} < \code{4/3}
        }
    }
}
\details{
  Also found in the SASmixed package
}
\source{
    Littel, R. C., Milliken, G. A., Stroup, W. W., and Wolfinger,
    R. D. (1996), \emph{SAS System for Mixed Models}, SAS Institute
    (Data Set 2.2(b)).
                                                                                
}
\keyword{datasets}
