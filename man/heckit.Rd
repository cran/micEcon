\name{heckit}
\alias{heckit}
\alias{coef.heckit}

\title{2-step Heckman (heckit) estimation}

\description{
   \code{heckit} performs a 2-step Heckman (heckit) estimation
   that corrects for non-random sample selection.
}

\usage{
   heckit( selection, formula, data, inst = NULL, print.level = 0 )
}

\arguments{
  \item{selection}{formula for the probit estimation (1st step) (see details).}
  \item{formula}{formula to be estimated (2nd step).}
  \item{data}{a data frame containing the data.}
  \item{inst}{an optional one-sided formula specifying instrumental variables
    for a 2SLS/IV estimation on the 2nd step.}
  \item{print.level}{this argument determines the level of printing which is
    done during the estimation process.  The default value of
    '0' means that no printing occurs, a value of '1' means that
    \code{heckit} reports what is currently done.}
}

\details{
   The endogenous variable of the argument 'selection' must have exactly
   two levels (e.g. 'FALSE' and 'TRUE', or '0' and '1').
   Only those observations are included in the second step estimation
   (argument 'formula'), where this variable equals the second element
   of its levels.
   By default the levels are sorted in increasing order
   ('FALSE' is before 'TRUE', and '0' is before '1').
   Hence, the second step is generally estimated for those
   observations, where the endogenous variable of the first step
   estimation equals 'TRUE' or '1'.
}

\value{
   \code{heckit} returns an object of class 'heckit' containing
   following elements:

   \item{coef}{estimated coefficients, standard errors, t-values and
   P-values.}

   \item{vcov}{variance covariance matrix of the estimated coefficients.}

   \item{probit}{object of class 'glm' that contains the results of the
   1st step (probit estimation).}

   \item{lm}{object of class 'lm' that contains the results of the
   2nd step (linear estimation). Note: the standard errors of this
   estimation are biased, because they do not account for the
   estimation of \eqn{\gamma} in the 1st step estimation
   (the correct standard errors are returned in \code{coef}.}

   \item{sigma}{the estimated \eqn{\sigma}, the standard error of the residuals.}

   \item{rho}{the estimated \eqn{\rho}, see Greene (2003, p. 784).}

   \item{invMillsRatio}{the inverse Mills Ratios calculated from the results of
   the 1st step probit estimation.}

   \item{imrDelta}{the \eqn{\delta}s calculated from the inverse Mills Ratios and
   the results of the 1st step probit estimation.}
}

\references{

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Johnston, J. and J. DiNardo (1997)
  \emph{Econometric Methods, Fourth Edition}, McGraw-Hill.

  Lee, L., G. Maddala and R. Trost (1980)
  Asymetric covariance matrices of two-stage probit and two-stage tobit
  methods for simultaneous equations models with selectivity.
  \emph{Econometrica}, 48, p. 491-503.

  Wooldridge, J. M. (2003)
  \emph{Introductory Econometrics: A Modern Approach, 2e}, Thomson South-Western.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{lm}}, \code{\link{glm}}, \code{\link{binomial}}}


\examples{
## Greene( 2003 ): example 22.8, page 786
data( Mroz87 )
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
greene <- heckit( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
summary( greene )        # print summary
summary( greene$probit ) # summary of the 1st step probit estimation
                         # this is Example 21.4, p. 681f
greene$sigma             # estimated sigma
greene$rho               # estimated rho

## Wooldridge( 2003 ): example 17.5, page 590
data( Mroz87 )
wooldridge <- heckit( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87 )
summary( wooldridge )        # summary of the 1st step probit estimation
                             # (Example 17.1, p. 562f) and 2nd step OLS regression
wooldridge$sigma             # estimated sigma
wooldridge$rho               # estimated rho

## example using random numbers
library( MASS )
nObs <- 1000
sigma <- matrix( c( 1, -0.7, -0.7, 1 ), ncol = 2 )
errorTerms <- mvrnorm( nObs, c( 0, 0 ), sigma )
myData <- data.frame( no = c( 1:nObs ), x1 = rnorm( nObs ), x2 = rnorm( nObs ),
   u1 = errorTerms[ , 1 ], u2 =  errorTerms[ , 2 ] )
myData$y <- 2 + myData$x1 + myData$u1
myData$s <- ( 2 * myData$x1 + myData$x2 + myData$u2 - 0.2 ) > 0
myData$y[ !myData$s ] <- NA
myOls <- lm( y ~ x1, data = myData)
summary( myOls )
myHeckit <- heckit( s ~ x1 + x2, y ~ x1, myData, print.level = 1 )
summary( myHeckit )

## example using random numbers with IV/2SLS estimation
library( MASS )
nObs <- 1000
sigma <- matrix( c( 1, 0.5, 0.1, 0.5, 1, -0.3, 0.1, -0.3, 1 ), ncol = 3 )
errorTerms <- mvrnorm( nObs, c( 0, 0, 0 ), sigma )
myData <- data.frame( no = c( 1:nObs ), x1 = rnorm( nObs ), x2 = rnorm( nObs ),
   u1 = errorTerms[ , 1 ], u2 = errorTerms[ , 2 ], u3 = errorTerms[ , 3 ] )
myData$w <- 1 + myData$x1 + myData$u1
myData$y <- 2 + myData$w + myData$u2
myData$s <- ( 2 * myData$x1 + myData$x2 + myData$u3 - 0.2 ) > 0
myData$y[ !myData$s ] <- NA
myHeckit <- heckit( s ~ x1 + x2, y ~ w, data = myData )
summary( myHeckit )  # biased!
myHeckitIv <- heckit( s ~ x1 + x2, y ~ w, data = myData, inst = ~ x1 )
summary( myHeckitIv ) # unbiased
}

\keyword{models}
