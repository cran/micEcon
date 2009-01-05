\name{coef.quadFuncEst}
\alias{coef.quadFuncEst}
\alias{vcov.quadFuncEst}

\title{Coefficients of a Quadratic Function}

\description{
   These methods return the coefficients and their covariance matrix
   from an estimated quadratic function.
}

\usage{
   \method{coef}{quadFuncEst}( object, \dots )

   \method{vcov}{quadFuncEst}( object, \dots )
}

\arguments{
  \item{object}{an object of class \code{quadFuncEst}.}
  \item{\dots}{currently ignored.}
}

\value{
The \code{coef} method returns a vector containing all (linearly independent)
coefficients of a quadratic function.

The \code{vcov} method returns the covariance matrix
of all (linearly independent) coefficients of a quadratic function.
}

\author{Arne Henningsen, \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{quadFuncEst}}}

\keyword{methods}
