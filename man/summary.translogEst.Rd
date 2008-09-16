\name{summary.translogEst}
\alias{summary.translogEst}
\alias{print.summary.translogEst}
\title{Summarizing the Estimation of a Translog Function}

\description{
   \code{summary.translogEst} summarizes the estimation results
   of a Translog Function.
}

\usage{
\method{summary}{translogEst}( object, ... )

\method{print}{summary.translogEst}( x, ... )
}

\arguments{
   \item{object}{an object of class \code{translogEst}.}
   \item{x}{an object of class \code{summary.translogEst}.}
   \item{...}{currently ignored.}
}

\value{
   \code{summary.translogEst} returns a list of class \code{summary.translogEst}
   that is currently the provided \code{object},
   but an element \code{coefTable} has been added
   and the class has been changed.
}

\seealso{\code{\link{translogEst}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\keyword{models}
