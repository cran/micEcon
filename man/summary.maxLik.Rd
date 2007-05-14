\name{summary.maxLik}
\alias{summary.maxLik}
\title{summary the Maximum-Likelihood estimation}
\description{
  Summary the Maximum-Likelihood estimation including standard errors
  and t-values.
}
\usage{
\method{summary}{maxLik}(object, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{object of class 'maxLik', most probably a result from
    Maximum-Likelihood estimation.}
  \item{\ldots}{currently not used.}
}
\value{
  object of class 'summary.maxLik' with following components:
  \item{type}{type of maximisation.}
  \item{iterations}{number of iterations.}
  \item{code}{code of success.}
  \item{message}{a short message describing the code.}
  \item{loglik}{the loglik value in the maximum.}
  \item{estimate}{numeric matrix, the first column contains the parameter
    estimates, the second the standard errors, third t-values and fourth
    corresponding probabilities.}
  \item{activePar}{logical vector, which parameters are treated as free.}
  \item{NActivePar}{number of free parameters.}
}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{maxLik}}}
\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, print.level=2)
summary(a)
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
summary(a)
}
\keyword{models}
