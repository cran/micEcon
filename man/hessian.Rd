\name{hessian}
\alias{hessian}
\alias{hessian.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Hessian matrix}
\description{
  This function extracts the Hessian of the M-estimator of 
  statistical model.  It should be supplied by the underlying
  optimisation algorithm, possibly using approximations.}
}
\usage{
hessian(x, \dots)
hessian.default(x, \dots)
}

\arguments{
  \item{x}{a M-estimator based statistical model}
  \item{\dots}{other arguments for methods}
}
\value{
  A numeric matrix, the Hessian of the model at the estimated parameter
  values.  If the maximum is flat, the Hessian
  is singular.  In that case you may want to invert only the
  non-singular part of the matrix.  You may also want to fix certain
  parameters (see \code{\link{activePar}}).
}
\author{Ott Toomet, \email{otoomet@econ.au.dk}}
\seealso{\code{\link{maxLik}}, \code{\link{activePar}}}
\examples{
# log-likelihood for normal density
# a[1] - mean
# a[2] - standard deviation
ll <- function(a) sum(-log(a[2]) - (x - a[1])^2/(2*a[2]^2))
x <- rnorm(1000) # sample from standard normal
ml <- maxLik(ll, start=c(1,1))
# ignore eventual warnings "NaNs produced in: log(x)"
summary(ml) # result should be close to c(0,1)
hessian(ml) # How the Hessian looks like
sqrt(-solve(hessian(ml))) # Note: standard deviations are on the diagonal
#
# Now run the same example while fixing a[2] = 1
mlf <- maxLik(ll, start=c(1,1), activePar=c(TRUE, FALSE))
summary(mlf) # first parameter close to 0, the second exactly 1.0
hessian(mlf)
# Now look at the Hessian.  Note that the second component has is not
# related to the curvature of the object function:
sqrt(-solve(hessian(mlf))) # works in current case, although you may get
                           # warnings.  
# now invert only the free parameter part of the Hessian
sqrt(-solve(hessian(mlf)[activePar(mlf), activePar(mlf)]))
# gives the standard deviation for the mean
}
\keyword{methods}
\keyword{optimize}
