\name{vcov-methods}
\alias{vcov.probit}
\alias{vcov.maxLik}

\title{'vcov' methods in Package 'micEcon'}
\description{
   Extract variance-covariance matrices.
}
\usage{
   \method{vcov}{probit}( object, ... )
   \method{vcov}{maxLik}( object, ... )
}
\arguments{
   \item{object}{an object of class \code{probit} or \code{maxLik}.}
   \item{\dots}{further arguments (currently ignored).}
}
\value{
   the estimated variance covariance matrix of the coefficients.
}
\author{
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet \email{otoomet@ut.ee}
}

\seealso{\code{\link[stat]{vcov}}, \code{\link{maxLik}}.}
\keyword{methods}
