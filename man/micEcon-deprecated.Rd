\name{micEcon-deprecated}
\alias{micEcon-deprecated}
\alias{tobit2}
\title{Deprecated Functions in micEcon package}
\description{
  These functions are provided for compatibility with older versions of
  micEcon only, and may be defunct as soon as the next release.
}
\usage{
tobit2(selection, formula, data = sys.frame(sys.parent()),
method="ml",
       start=NULL, print.level=0,
       y1=FALSE, z=FALSE, y2=FALSE, x=FALSE, model=FALSE,
       ...)
}
\arguments{
  \item{selection}{a symbolic formula for the selection equation.  The
    response must be binary variable where 0 corresponds to invisibility
    of the outcome and 1 to visibility.}
  \item{formula}{a symbolic formula for the equation of interest.}
  \item{data}{an optional data frame containing the variables in
    the model. By default the variables are
    taken from environment(formula), typically the environment from
    which the function is called}
  \item{method}{the method of calculation: Maximum Likelihood ("\code{ml}") or
    Heckman 2-step (\code{"2step"})}
  \item{start}{numeric, initial values of the parameters.  The order is
    as follows: gamma (selection equation), beta (the equation of
    interest), sigma, rho (distribution of the error terms) for tobit2;}
  \item{print.level}{information about calculations.  0 -- nothing is
    printed, bigger numbers give more information.
    \code{print.level} is sent further to the maximisation algorithm,
    see \code{\link{maxNR}}.}
  \item{y1, z, y2, x, model}{logicals. If TRUE the corresponding
    components of the fit (the selection response, the selection model
    matrix, the equation response, the equation model matrix and both
    model frames) are returned.}
  \item{\dots}{further arguments to the maximisation algorithm,
    see \code{\link{maxNR}}.}
}
\details{
  The original help page for these functions is often
  available at \code{help("oldName-deprecated")} (note the quotes).

  Function \code{tobit2} are for historical
  compatibility with old versions of micEcon.  Use
  \code{\link{selection}} instead.
}
\seealso{
  \code{\link{Deprecated}}, \code{\link{selection}}
}
\keyword{misc}
