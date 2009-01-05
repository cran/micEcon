\name{coef.translogEst}
\alias{coef.translogEst}
\alias{vcov.translogEst}

\title{Coefficients of a Translog Function}

\description{
   These methods return the coefficients and their covariance matrix
   from an estimated translog function.
}

\usage{
   \method{coef}{translogEst}( object, \dots )

   \method{vcov}{translogEst}( object, \dots )
}

\arguments{
  \item{object}{an object of class \code{translogEst}.}
  \item{\dots}{currently ignored.}
}

\value{
The \code{coef} method returns a vector containing all (linearly independent)
coefficients of a translog function.

The \code{vcov} method returns the covariance matrix
of all (linearly independent) coefficients of a translog function.
}

\author{Arne Henningsen, \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{translogEst}}}

\keyword{methods}
