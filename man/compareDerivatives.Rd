\name{compareDerivatives}
\alias{compareDerivatives}
\title{function to compare analytic and numeric derivatives}
\description{
  This function compares analytic and numerical derivative and prints a few
  diagnostics.  It is intended for testing pre-programmed derivative
  routines
  for maximisation algorithms.
}
\usage{
compareDerivatives(f, grad, hess=NULL, t0, eps=1e-6, ...)
}
\arguments{
  \item{f}{function to be differentiated. The parameter (vector) of interest
    must be the first argument. The function may return a vector.}
  \item{grad}{function returning the analytic gradient.
    Must use the same set of parameters as \code{f}.
    If \code{f} is a vector-valued function, grad must return a matrix where the
    number of rows equals the number of components of \code{f}, and the number
    of columns must equal to the number of components in \code{t0}.}
  \item{hess}{function returning the analytic hessian.
    If present, hessian matrices are
    compared too.  Only appropriate for scalar-valued functions.}
  \item{t0}{parameter vector indicating the point at which the derivatives
    are compared. The derivative is taken with respect to this
    vector.}
  \item{eps}{numeric. Step size for numeric differentiation. Central
    derivative is used.}
  \item{\dots}{further arguments to \code{f}, \code{grad} and \code{hess}.}
}
\details{
  For every component of \code{f}, the parameter value, analytic and numeric
  derivative and their relative difference (analytic - numeric)/analytic
  are printed.  If analytic derivatives are correct and the function is
  sufficiently smooth, expect the relative differences to be less than 1e-7.
}
\author{Ott Toomet \email{siim@obs.ee}}
\seealso{\code{\link{numericGradient}}}

\examples{
## A simple example with sin(x)' = cos(x)
f <- sin
compareDerivatives(f, cos, t0=1)
##
## Example of log-likelihood of normal density.  Two-parameter
## function.
x <- rnorm(100, 1, 2) # generate rnorm x
l <- function(b) sum(log(dnorm((x-b[1])/b[2])/b[2]))
              # b[1] - mu, b[2] - sigma
gradl <- function(b) {
c(sum(x - b[1])/b[2]^2,
sum((x - b[1])^2/b[2]^3 - 1/b[2]))
}
compareDerivatives(l, gradl, t0=c(1,2))
}
\keyword{math}
\keyword{utilities}
