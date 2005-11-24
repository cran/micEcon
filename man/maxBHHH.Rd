\name{maxBHHH}
\alias{maxBHHH}
\title{BHHH maximisation}
\description{
This function is essentially \code{\link{maxNR}} where the hessian is
approximated by the information equation.}
}
\usage{
maxBHHH(fn, grad = NULL, hess=NULL, theta, print.level = 0, iterlim = 100,
   ...)
}

\arguments{
  \item{fn}{function to be maximised, must return vector of
    observation-specific likelihood values.  Those are summed for a
    single likelihood value.  If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.}
  \item{grad}{gradient of the function.  If \code{NULL}, numeric
    gradient is used.  It must return a matrix, where
    rows corresponds to the gradient vectors of individual
    observations.  Note that this corresponds to
    \code{t(numericGradient(fn))}, not \code{numericGradient(fn)}.
    It is summed over
    observations in order to get a single gradient vector.}
  \item{hess}{Hessian of the function.  Not used, for argument
    compatibility reasons.}
  \item{theta}{initial values for the parameters to be optimized over.}
  \item{print.level}{a larger number prints more working information.}
  \item{iterlim}{maximum number of iterations.}
  \item{\dots}{further arguments for \code{\link{maxNR}}, \code{fn} and \code{grad}.}
}
\details{
  \code{maxBHHH} uses information equality in order to approximate the
  hessian of the log-likelihood function.  This must be average (summed)
  over individual observations.  Hence we need log-likelihood and
  gradient by observations, hessian is calculated as \code{-t(gradient)
    \%*\% gradient}.
}
\value{
  Object of class "maximisation":
  \item{type}{Character string: "BHHH maximisation"}
  Plus the other components inherited from \code{\link{maxNR}}.
}

\author{Ott Toomet \email{otoomet@ut.ee}}

\section{Warning }{As \code{maxBHHH} uses likelihood-specific
  information equation, it is only suitable for maximising log-likelihood!}

\seealso{\code{\link{maxNR}}, \code{\link{maxBFGS}},
  \code{\link{optim}}, \code{\link{nlm}}.}

\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
## Estimate with numeric gradient and hessian
a <- maxBHHH(loglik, theta=1, print.level=2)
summary(a)
## Estimate with analytic gradient
a <- maxBHHH(loglik, gradlik, theta=1)
summary(a)
}
\keyword{optimize}
