\name{returnMessage}
\alias{returnMessage}
\alias{returnMessage.default}
\alias{returnMessage.maximisation}
\alias{returnMessage.maxLik}
\title{Information about the optimisation process}
\description{
  This function returns a short message, summarising the outcome of the
  statistical process, typically optimisation.  The message should describe either
  the type of the convergence,
  or the problem.  \code{returnMessage} is a generic function, with methods
  for various optimisation algorithms.
}
\usage{
returnMessage(x, ...)
returnMessage.maximisation(x, ...)
returnMessage.maxLik(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object, should orginate from an optimisation problem}
  \item{...}{further arguments to other methods}.
}
\details{
  The default methods returns component \code{returnMessage}.
}
\value{
  Character string, the message describing the success or failure of the
  statistical procedure.
}
\author{Ott Toomet, \email{otoomet@ut.ee}}
\seealso{\code{\link{returnCode}}, \code{\link{maxNR}}}
\examples{
## maximise the exponential bell
f1 <- function(x) exp(-x^2)
a <- maxNR(f1, start=2)
returnCode(a) # should be success (1 or 2)
## Now try to maximise log() function
f2 <- function(x) log(x)
a <- maxNR(f2, start=2)
returnCode(a) # should give a failure (4)
}
\keyword{methods}
\keyword{utilities}
