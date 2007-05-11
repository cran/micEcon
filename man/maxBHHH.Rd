\name{maxBHHH}
\alias{maxBHHH}
\title{BHHH maximisation}
\description{
This function is essentially \code{\link{maxNR}} where the hessian is
approximated by the outer product of the gradient vectro.  Only suitable
for likelihood maximisation.}
}
\usage{
maxBHHH(fn, grad = NULL, hess=NULL, start, print.level = 0, iterlim = 100,
   ...)
}

\arguments{
  \item{fn}{log-likelihood function to be maximised, must have the
    parameter vector as
    the first argument and must return vector
    of
    observation-specific likelihood values.  Those are summed for a
    single likelihood value.  If the parameters are out of range, \code{fn} should
    return \code{NA}.}
  \item{grad}{gradient of the log-likelihood function.  Must have the parameter vector as
    the first argument.  If \code{NULL}, numeric
    gradient is used.  It must return a matrix, where
    rows corresponds to the gradient vectors of individual
    observations.  This matrix is summed over
    observations in order to get a single gradient vector.  It's matrix
    product is used for approximation of the Hessian.}
  \item{hess}{Hessian of the function.  Not used by \code{maxBHHH}, for
    compatibility with \code{\link{maxNR}}.}
  \item{start}{initial values for the parameters to be optimized over.}
  \item{print.level}{a larger number prints more working information.}
  \item{iterlim}{maximum number of iterations.}
  \item{\dots}{further arguments for \code{\link{maxNR}}, \code{fn} and \code{grad}.}
}
\details{
  \code{maxBHHH} uses information equality in order to approximate the
  Hessian of the log-likelihood function.  Hence we have to calculate
  log-likelihood
  and
  gradient by individual observations.  Hessian is approximated as sum
  of outer products of the gradients of individual observations, or, in
  the matrix form, \code{-t(gradient)
    \%*\% gradient}.
}
\value{
  Object of class "maximisation":
  \item{type}{Character string: "BHHH maximisation"}
  Plus the other components inherited from \code{\link{maxNR}}.
}

\author{Ott Toomet \email{otoomet@ut.ee}}

\section{Warning}{As \code{maxBHHH} uses likelihood-specific
  information equation, it is only suitable for maximising log-likelihood!}

\seealso{\code{\link{maxNR}}, \code{\link{maxBFGS}},
  \code{\link{optim}}, \code{\link{nlm}}.}

\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
## Estimate with numeric gradient and hessian
a <- maxBHHH(loglik, start=1, print.level=2)
summary(a)
## Estimate with analytic gradient
a <- maxBHHH(loglik, gradlik, start=1)
summary(a)
}
\keyword{optimize}
