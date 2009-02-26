\name{elas}
\alias{elas}
\alias{elasticities}
\alias{elas.default}

\title{Calculating and returning elasticities}

\description{
   These functions calculate and return
   elasticities of microeconomic models.
   \code{elasticities} is an alias for \code{elas}.
}

\usage{
   elas( object, \dots )
   elasticities( object, \dots )
   \method{elas}{default}( object, \dots )
}

\arguments{
  \item{object}{a microeconomic model}
  \item{\dots}{further arguments for methods}
}

\details{
   This is a generic function.
   The default method just returns the element
   \code{elas} from \code{object}.
}

\author{Arne Henningsen}

\seealso{\code{\link{elas.aidsEst}}}

\keyword{methods}
