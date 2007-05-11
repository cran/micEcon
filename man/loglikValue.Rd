\name{loglikValue}
\alias{loglikValue}
\alias{loglikValue.default}
\alias{loglikValue.maxLik}
\title{Return the log likelihood value}
\description{
  Return the log likelihood value.  Objects should be models,
  estimated by Maksimum Likelihood, or other statistical models where
  the likelihood makes sense.
}
\usage{
loglikValue(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object, usually a model estimated with Maximum
    Likelihood}
  \item{...}{additional arguments to methods}
}
\details{
  The default method returns the component \code{loglik} 
}
\value{
  A single numeric, log likelihood of the estimated model
}
\author{Ott Toomet \email{otoomet@ut.ee}}
\seealso{\code{\link{maxLik}}}
\examples{
## We estimate a tobit-2 model and print the loglik value
library(mvtnorm)
N <- 500
   vc <- diag(2)
   vc[2,1] <- vc[1,2] <- -0.7
   eps <- rmvnorm(N, rep(0, 2), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo <- runif(N)
   yo <- (xo + eps[,2])*(ys > 0)
   a <- selection(ys~xs, yo ~xo)
loglikValue(a) # this values is also printed by 'summary'
}
\keyword{methods}
