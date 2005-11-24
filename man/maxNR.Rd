\name{maxNR}
\alias{maxNR}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Newton-Raphson maximisation}
\description{
  unconstrained maximisation algorithm based on Newton-Raphson method.
}
\usage{
maxNR(fn, grad = NULL, hess = NULL, theta, print.level = 0,
      tol = 1e-06, gradtol = 1e-06, steptol = 1e-06, lambdatol = 1e-06,
      qrtol = 1e-10,
      iterlim = 15,
      constPar = NULL, activePar=rep(TRUE, NParam), ...)
}
\arguments{
  \item{fn}{function to be maximised.  In order to use numeric gradient
    and BHHH method, \code{fn} must return vector of
    observation-specific likelihood values.  Those are summed by maxNR
    if necessary.  If the parameters are out of range, \code{fn} should
    return \code{NA}.  See details for constant parameters.}
  \item{grad}{gradient of the function.  If \code{NULL}, numeric
    gradient is used.  It must return a gradient vector, or matrix where
    \emph{columns} correspond to individual parameters.  Note that this corresponds to
    \code{t(numericGradient(fn))}, not \code{numericGradient(fn)}.  For
    BHHH method,
    the rows must correspond to the likelihood gradients of the individual
    observations, it will be summed over
    observations in order to get a single gradient vector.  Similar
    observation-wise matrix is allowed for \code{maxNR} too, it will be
    simply summed over observations.}
  \item{hess}{hessian matrix of the function.  If missing, numeric
    hessian, based on gradient, is used.  }
  \item{theta}{initial value for the parameter vector}
  \item{print.level}{this argument determines the level of
    printing which is done during the minimization process. The default
    value of 0 means that no printing occurs, a value of 1 means that
    initial and final details are printed and a value of 2 means that
    full tracing information for every iteration is printed.
  }
  \item{tol}{stopping condition.  Stop if the absolute difference
    between successive iterations less than \code{tol}, return
    \code{code=2}.}
  \item{gradtol}{stopping condition.  Stop if the norm of the gradient
    less than \code{gradtol}, return \code{code=1}.}
  \item{steptol}{stopping/error condition.  If not find
    acceptable/higher value of the function with \code{step=1}, step is
    divided by 2 and tried again.  This is repeated until \code{step <
      steptol}, then \code{code=3} is returned.}
  \item{lambdatol}{control whether the hessian is treated as negative
    definite.  If the
    largest of the eigenvalues of the hessian is larger than
    \code{-lambdatol}, a suitable diagonal matrix is subtracted from the
    hessian (quadratic hill-climbing).}
  \item{qrtol}{?}
  \item{iterlim}{stopping condition.  Stop if more than \code{iterlim}
    iterations, return \code{code=4}.}
  \item{constPar}{index vector.  Which of the parameters must be treated
    as constants.}
  \item{activePar}{logical vector, which parameters are treated as free
    (resp constant)}
  \item{\dots}{further argument to \code{fn}, \code{grad} and
    \code{hess}.}
}
\details{
The algorithm can treat constant parameters and related changes of
parameter values.  Constant parameters are useful if a parameter value
is converging toward the edge of support or for testing.  

One way is to put
\code{constPar} to non-NULL.  Second possibility is to signal by
\code{fn} which parameters are constant and change the values of the
parameter vector.  The value of \code{fn} may have following attributes:
\item{constPar}{index vector.  Which parameters are redefined to
  constant}
\item{constVal}{numeric vector.  Values of the constant parameters.}
\item{newVal}{a list with following components:
  \item{index}{which parameters will have a new value}
  \item{val}{the new value of parameters}
}
The difference between \code{constVal} and \code{newVal} is that the
latter parameters are not set to constants.  If the attribute
\code{newVal} is present, the new function value is allowed to be below
the previous one.
}
\value{
  list of class "maximisation" which following components:
  \item{maximum}{\code{fn} value at maximum (the last calculated value
    if not converged).}
  \item{estimate}{estimated parameter value.}
  \item{gradient}{last gradient value which was calculated.  Should be
    close to 0 if normal convergence.}
  \item{hessian}{hessian at the maximum (the last calculated value if
    not converged).  May be used as basis for variance-covariance
    matrix.}
  \item{code}{return code:
    \item{1}{gradient close to zero (normal convergence).}
    \item{2}{successive function values within tolerance limit (normal
      convergence).}
    \item{3}{last step could not find higher value (probably not
      converged).}
    \item{4}{iteration limit exceeded.}
    \item{100}{Initial value out of range.}
  }
  \item{message}{a short message, describing \code{code}.}
  \item{last.step}{list describing the last unsuccessful step if
    \code{code=3} with following components:
    \item{theta0}{previous parameter value}
    \item{f0}{\code{fn} value at \code{theta0}}
    \item{theta1}{new parameter value.  Probably close to the previous
      one.}
  }
  \item{activePar}{logical vector, which parameters are not constants.}
  \item{iterations}{number of iterations.}
  \item{type}{character string, type of maximisation.}
}
\references{W. Greene: "Advanced Econometrics"; S.M. Goldfeld,
  R.E. Quandt: "Nonlinear Methods in Econometrics".  Amsterdam,
  North-Holland 1972.}
\author{Ott Toomet \email{siim@obs.ee}}
\note{Newton-Raphson algorithm with analytic gradient and hessian
  supplied should converge with no more than 20 iterations on a
  well-specified model.}
\seealso{\code{\link{nlm}} for Newton-Raphson optimisation,
  \code{\link{optim}} for different gradient-based optimisation
  methods.}
\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) sum(log(theta) - theta*t)
## Note the log-likelihood and gradient are summed over observations
gradlik <- function(theta) sum(1/theta - t)
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxNR(loglik, theta=1, print.level=2)
summary(a)
## You would probably prefer 1/mean(t) instead ;-)
## Estimate with analytic gradient and hessian
a <- maxNR(loglik, gradlik, hesslik, theta=1)
summary(a)
}
\keyword{optimize}
