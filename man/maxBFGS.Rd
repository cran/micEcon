\name{maxBFGS}
\alias{maxBFGS}
\alias{maxSANN}
\alias{maxNM}

\title{BFGS, SANN and Nelder-Mead Maximization}
\description{
  These functions are wrappers for \code{\link{optim}} where the arguments are
  compatible with \code{\link{maxNR}}
}
\usage{
maxBFGS(fn, grad = NULL, hess=NULL, start, print.level = 0, iterlim = 200,
   tol = 1e-08, reltol=tol, ... )
maxSANN(fn, grad = NULL, hess = NULL, start, print.level = 0, iterlim =
   10000, tol = 1e-08, reltol=tol, temp = 10, tmax = 10, parscale = rep(1, length = length(start)), ...)
maxNM(fn, grad = NULL, hess = NULL, start, print.level = 0, iterlim =
   500, tol = 1e-08, reltol=tol, parscale = rep(1, length = length(start)), alpha = 1, beta = 0.5, gamma = 2, ...)
}

\arguments{
  \item{fn}{function to be maximised.  Must have the parameter vector as
    the first argument.  In order to use numeric gradient
    and BHHH method, \code{fn} must return vector of
    observation-specific likelihood values.  Those are summed by maxNR
    if necessary.  If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.}
  \item{grad}{gradient of the function.  Must have the parameter vector as
    the first argument.  If \code{NULL}, numeric
    gradient is used (only maxBFGS uses gradient).  Gradient may return
    a matrix, where columns correspond to the parameters and rows to the
    observations (useful for maxBHHH).  The columns are summed internally.}
  \item{hess}{Hessian of the function.  Not used by any of these methods, for
    compatibility with \code{\link{maxNR}}.}
  \item{start}{initial values for the parameters.}
  \item{print.level}{a larger number prints more working information.}
  \item{iterlim}{maximum number of iterations.}
  \item{tol, reltol}{the relative convergence tolerance (see
    \code{\link{optim}}).  \code{tol} is for compatibility with \code{\link{maxNR}}.}
  \item{temp}{controls the '"SANN"' method. It is the starting
          temperature for the cooling schedule. Defaults to '10'.}
  \item{tmax}{is the number of function evaluations at each temperature
          for the '"SANN"' method. Defaults to '10'. (see
    \code{\link{optim}})}
  \item{parscale}{A vector of scaling values for the parameters.
          Optimization is performed on 'par/parscale' and these should
          be comparable in the sense that a unit change in any element
          produces about a unit change in the scaled value. (see
    \code{\link{optim}})}
  \item{alpha, beta, gamma}{Scaling parameters for the
          '"Nelder-Mead"' method. 'alpha' is the reflection factor
          (default 1.0), 'beta' the contraction factor (0.5) and
          'gamma' the expansion factor (2.0). (see
    \code{\link{optim}})}
  \item{\dots}{further arguments for \code{fn} and \code{grad}.}
}
\value{
  Object of class "maximisation":
  \item{maximum}{value of \code{fn} at maximum.}
  \item{estimate}{best set of parameters found.}
  \item{gradient}{gradient at parameter value \code{estimate}.}
  \item{hessian}{value of Hessian at optimum.}
  \item{code}{integer. Success code, 0 is success (see
    \code{\link{optim}}).}
  \item{message}{character string giving any additional information
    returned by the optimizer, or NULL.}
  \item{iterations}{two-element integer vector giving the number of
    calls to \code{fn} and \code{gr}, respectively.
    This excludes those calls needed to
    compute the Hessian, if requested, and any calls to \code{fn} to compute a
    finite-difference approximation to the gradient.}
  \item{type}{character string "BFGS maximisation".}
}
\author{Ott Toomet \email{otoomet@ut.ee}}

\seealso{\code{\link{optim}}, \code{\link{nlm}}, \code{\link{maxNR}}, \code{\link{maxBHHH}}.}

\examples{
# Maximum Likelihood estimation of the parameter of Poissonian distribution
n <- rpois(100, 3)
loglik <- function(l) n*log(l) - l - lfactorial(n)
# we use numeric gradient
summary(maxBFGS(loglik, start=1))
# you would probably prefer mean(n) instead of that ;-)
# Note also that maxLik is better suited for Maximum Likelihood
}
\keyword{optimize}
