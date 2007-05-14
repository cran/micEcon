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
   elas( object, ... )
   elasticities( object, ... )
   \method{elas}{default}( object, ... )
}

\arguments{
  \item{object}{a microeconomic model}
  \item{...}{further arguments for methods}
}

\details{
   This is a generic function.
   The default method just returns the element
   \code{elas} from \code{object}.
}

\author{Arne Henningsen, \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{elas.aidsEst}}}

\keyword{methods}
