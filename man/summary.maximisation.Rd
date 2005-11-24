\name{summary.maximisation}
\alias{summary.maximisation}
\title{Summary method for maximisation/minimisation}
\description{
  Selects the appropriate information from maximisation algorithm,
  calculates standard errors and
  \eqn{t}-values, and presents in a suitable form.
}
\usage{
summary.maximisation( object, hessian=FALSE, unsucc.step=FALSE, ... )
}
\arguments{
  \item{object}{optimisation result, object of class
    \code{maximisation}.  See \code{\link{maxNR}}.}
  \item{hessian}{logical, whether to display Hessian matrix.}
  \item{unsucc.step}{logical, whether to describe last unsuccesful step
    if \code{code} == 3}
  \item{\ldots}{currently not used.}
}
\value{
  Object of class \code{summary.maximisation}, intended to print with
  corresponding print method.  There are following components:
  \item{type}{type of maximisation.}
  \item{iterations}{number of iterations.}
  \item{code}{exit code (see \code{\link{maxNR}}.)}
  \item{message}{a brief message, explaining code.}
  \item{unsucc.step}{description of last unsuccessful step, only if
    requested and \code{code} == 3}
  \item{maximum}{function value at maximum}
  \item{estimate}{matrix with following columns:
    \itemize{
      \item{results:}{coefficient estimates at maximum}
      \item{gradient:}{estimated gradient at maximum}
    }
  }
  \item{hessian}{estimated hessian at maximum, only if requested}
}
\author{Ott Toomet \email{siim@obs.ee}}
\section{WARNING }{Probabilities of t-values depend on several
  assumptions which may not be appropriate in the particular case!}
\seealso{\code{\link{maxNR}}}
\examples{
## minimize a 2D quadratic function:
f <- function(b) { x <- b[1]; y <- b[2];
    val <- (x - 2)^2 + (y - 3)^2
    attr(val, "gradient") <- c(2*x - 4, 2*y - 6)
    attr(val, "hessian") <- matrix(c(2, 0, 0, 2), 2, 2)
    val
}
## Note that NR finds the minimum of a quadratic function with a single
## iteration.  Use c(0,0) as initial value.  Note that standard
## deviations, t-values and probabilities are not appropriate here.
result1 <- maxNR( f, theta = c(0,0) ) 
summary( result1 )
## Now use c(1000000, -777777) as initial value and ask for hessian
result2 <- maxNR( f, theta = c( 1000000,-777777 ) ) 
summary( result2 )
}
\keyword{methods}% at least one, from doc/KEYWORDS
\keyword{print}
