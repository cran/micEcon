\name{maximisationType}
\alias{maximisationType}
\alias{maximisationType.default}
\alias{maximisationType.maximisation}
\alias{maximisationType.MLEstimate}
\title{Type of minimisation/maximisation}
\description{
  Returns the type of optimisation.  It should be returned
  by the optimisation routine.
}
\usage{
  maximisationType(x)
}
\arguments{
  \item{x}{object of class 'maximisation' or another object which
    involves numerical optimisation.
  }
}
\value{
  A text message, describing the involved optimisation algorithm
}
\author{Ott Toomet, \email{otoomet@ut.ee}}
\seealso{\code{\link{maxNR}}}
\examples{
## maximise two-dimensional exponential hat.  Maximum is at c(2,1):
f <- function(a) exp(-(a[1] - 2)^2 - (a[2] - 1)^2)
m <- maxNR(f, start=c(0,0))
summary(m)
maximisationType(m)
## Now try BHHH maximisation.  It shouldn't work but does.
m <- maxBHHH(f, start=c(0,0))
summary(m)
maximisationType(m)
}
\keyword{optimize}
\keyword{methods}
