\name{maxLik}
\alias{maxLik}
\alias{coef.maxLik}

\title{Maximum likelihood estimation}
\description{
  This is just a wrapper for maximisation routines which return object
  of calss "maxLik".  Corresponding methods can correctly handle the
  likelihood-specific properties of the estimate, e.g. the fact that
  inverse of negative hessian is the variance-covariance matrix.
}
\usage{
maxLik(logLik, grad = NULL, hess = NULL, theta, method = "Newton-Raphson", ...)
}
\arguments{
  \item{logLik}{log-likelihood function.  Must have the parameter vector
    as the first argument.  Must return either a single log-likelihood
    value or a numeric vector where each component is log-likelihood
    corresponding to individual observations.}
  \item{grad}{gradient of log-likelihood.  Must have the parameter
    vector as the first argument.  Must return either single gradient
    vector with length equal to the number of parameters, or a matrix
    where each row corresponds to gradient vector of individual
    observations.  If \code{NULL}, numeric gradient will be used.}
  \item{hess}{hessian of log-likelihood.  Must have the parameter
    vector as the first argument.  Must return a square matrix.  If
    \code{NULL}, numeric gradient will be used.}
  \item{theta}{numeric vector, initial value of parameters.}
  \item{method}{maximisation method, currently only "Newton-Rapshon".}
  \item{\dots}{further arguments for the maximisation routine.}
}
\value{
  object of class 'maxLik' which inherits from class 'maximisation'.
  Components are identical to those of class 'maximisation',
  see \code{\link{maxNR}}.
}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{maxNR}}, \code{\link{nlm}} and \code{\link{optim}}
  for different non-linear optimisation routines.}
\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, theta=1, print.level=2)
summary(a)
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, theta=1)
summary(a)
}
\keyword{optimize}
