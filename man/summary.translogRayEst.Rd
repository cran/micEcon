\name{summary.translogRayEst}
\alias{summary.translogRayEst}
\alias{print.summary.translogRayEst}
\title{Summarizing the Estimation of a Translog Ray Function}

\description{
   \code{summary.translogEst} summarizes the estimation results
   of a translog ray function.
}

\usage{
\method{summary}{translogRayEst}( object, \dots )

\method{print}{summary.translogRayEst}( x, \dots )
}

\arguments{
   \item{object}{an object of class \code{translogRayEst}.}
   \item{x}{an object of class \code{summary.translogRayEst}.}
   \item{\dots}{currently ignored.}
}

\value{
   \code{summary.translogRayEst} returns
   a list of class \code{summary.translogRayEst}
   that is currently the provided \code{object},
   but an element \code{coefTable} has been added
   and the class has been changed.
}

\seealso{\code{\link{translogRayEst}}.}

\author{Arne Henningsen}

\keyword{models}
