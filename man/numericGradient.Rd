\name{numericGradient}
\alias{numericGradient}
\alias{numericHessian}
\alias{numericNHessian}

\title{Functions to Calculate Numeric Derivatives}
\description{
  Calculate (central) numeric gradient and Hessian.
  \code{numericGradient} accepts vector-valued functions.
}
\usage{
numericGradient(f, t0, eps=1e-06, \dots)
numericHessian(f, grad=NULL, t0, eps=1e-06, \dots)
numericNHessian(f, t0, eps=1e-6, \dots)

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{f}{function to be differentiated.  The first argument must be
    the parameter vector with respect to which it is differentiated.}
  \item{grad}{function, gradient of f}
  \item{t0}{vector, the value of parameters}
  \item{eps}{numeric, the step for numeric differentiation}
  \item{\dots}{furter arguments for \code{f}}
}
\details{
  \code{numericGradient} numerically differentiates a (vector valued)
  function with respect to it's (vector valued) argument.  If the
  functions value is a \code{NVal * 1} vector and the argument is
  \code{Npar * 1} vector, the resulting gradient is a \code{NVal * NPar}
  matrix. 

  \code{numericHessian} checks whether a gradient function is present
  and calculates a gradient of the gradient (if present), or full
  numeric Hessian (\code{numericNHessian}) if \code{grad} is \code{NULL}.
}
\value{
  Matrix.  For \code{numericGradient}, the number of rows is equal to the
  length of the function value vector, and the number of columns is
  equal to the length of the parameter vector.

  For the \code{numericHessian}, both numer of rows and columns is
  equal to the length of the parameter vector.
}
\section{Warning}{
  Be careful when using numerical differentiation in optimisation
  routines.  Although quite precise in simple cases, they may work very
  poorly in more complicated cases.
  }
\author{Ott Toomet \email{siim@obs.ee}}
\seealso{\code{\link{compareDerivatives}}, \code{\link{deriv}}}

\examples{
# A simple example with Gaussian bell
f0 <- function(t0) exp(-t0[1]^2 - t0[2]^2)
numericGradient(f0, c(1,2))
numericHessian(f0, t0=c(1,2))

# An example with the analytic gradient
gradf0 <- function(t0) -2*t0*f0(t0)
numericHessian(f0, gradf0, t0=c(1,2))
# The results should be similar as in the previous case

# The central numeric derivatives have usually quite a high precision
compareDerivatives(f0, gradf0, t0=1:2)
# The differenc is around 1e-10
}
\keyword{math}
\keyword{utilities}
