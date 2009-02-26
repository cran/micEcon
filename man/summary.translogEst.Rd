\name{summary.translogEst}
\alias{summary.translogEst}
\alias{print.summary.translogEst}
\title{Summarizing the Estimation of a Translog Function}

\description{
   \code{summary.translogEst} summarizes the estimation results
   of a Translog Function.
}

\usage{
\method{summary}{translogEst}( object, \dots )

\method{print}{summary.translogEst}( x, \dots )
}

\arguments{
   \item{object}{an object of class \code{translogEst}.}
   \item{x}{an object of class \code{summary.translogEst}.}
   \item{\dots}{currently ignored.}
}

\value{
   \code{summary.translogEst} returns a list of class \code{summary.translogEst}
   that is currently the provided \code{object},
   but an element \code{coefTable} has been added
   and the class has been changed.
}

\seealso{\code{\link{translogEst}}.}

\author{Arne Henningsen}

\keyword{models}
