\name{returnCode}
\alias{returnCode}
\alias{returnCode.default}
\title{Return code for optimisation and other objects}
\description{
  This function gives the return code of various optimisation objects.
  The return code gives a brief information about the success or
  problems, occured during the optimisation (see documentation for the
  corresponding function).
}
\usage{
returnCode(x, ...)
\method{returnCode}{default}(x, ...)
}
\arguments{
  \item{x}{object, usually an estimator, achieved by optimisation}
  \item{...}{further arguments for other methods}
}
\details{
  The default methods returns component \code{returnCode}.
}
\value{
  Integer, the success code of optimisation procedure.  However,
  different optimisation routines may define it in a different way.
}
\author{Ott Toomet, \email{otoomet@ut.ee}}
\seealso{\code{\link{returnMessage}}, \code{\link{maxNR}}}
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
