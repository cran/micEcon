\name{tobit2}
\alias{tobit2}
\title{function to calculate Heckman-type selection models}
\description{
  Heckman-type selection models are regressions with at least two latent
  equation and multivariate normal error term.
  Tobit-2 model is one particular type of those models (Amemiya 1985).
}
\usage{
tobit2(selection, formula, data = sys.frame(sys.parent()), method="ml",
       b0=NULL, print.level=0,
       y1=FALSE, z=FALSE, y2=FALSE, x=FALSE, model=FALSE,
       ...)
}
\arguments{
  \item{selection}{a symbolic formula for the selection equation}
  \item{formula}{a symbolic formula for the equation of interest}
  \item{data}{an optional data frame containing the variables in
    the model. By default the variables are
    taken from environment(formula), typically the environment from
    which the function is called}
  \item{method}{the method of calculation: Maximum Likelihood (\code{ml}) or
    Heckman 2-step (\code{2step})}
  \item{b0}{numeric, initial values of the parameters.  The order is
    as follows: gamma (selection equation), beta (the equation of
    interest), sigma, rho (distribution of the error terms)}
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
  Tobit 2 models have two latent variable regressions with jointly
  normal errors.  If selection variable \eqn{y_1 = 0}{y1 = 0}, the
  formula response \eqn{y_2}{y2} is not observable; if \eqn{y_1 = 1}{y1
    = 1}, then \eqn{y_2 = y_2^*}{y2 = y2*}.
  Values for unobservable selection are undefined, they are commonly
  assumed to be equal to 0.  In the current function the values must be
  there but they are not used.
}
\value{
  a list of class \code{tobit2}:
  \item{maximum}{likelihood value at maximum}
  \item{estimate}{(named) vector of estimates}
  \item{twoStep}{object of class "heckit", results of two-step
    estimation}
  \item{NParam}{total number of parameters}
  \item{NObs}{total number of observations}
  \item{N1}{number of observations where the response was not observed}
  \item{N2}{number of observations where the response was observed}
  \item{NZ}{number of variables in the selection equation}
  \item{NX}{number of variables in the main formula}
  \item{df}{resulting degrees of freedom}
  \item{call, terms1, terms2}{call and terms of the model}
  Plus the components of \code{\link{maxLik}}.
}
\references{Amemiya 1985}
\author{Ott Toomet \email{siim@obs.ee}}
\section{WARNING}{These models rely heavily on the parametric
  assumptions.  Do not use them unless the exclusion restriction holds!}
\seealso{\code{\link{heckit}} for the two-step estimation and
   \code{\link{probit}} for probit models}
\examples{
## uniform observable covariates (note that exclusion restriction is
## fulfilled)!
myData <- data.frame( z = runif(500), x = runif(500) )
## latent selection variable
myData$y1s <- 2 * myData$z + rnorm(500)
myData$y1  <- myData$y1s > 0
## latent regressions
myData$y2s <- myData$x + 0.5 * rnorm(500)
myData$y2  <- ifelse( myData$y1 == 1, myData$y2s, 0)
plot( myData$z, myData$y2 )
a <- tobit2( y1~z, y2~x, data = myData )
summary(a)
}
\keyword{models}
\keyword{regression}
