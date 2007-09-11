\name{selection}
\alias{selection}
\alias{heckit}
\alias{print.selection}
\title{Heckman-style selection models}
\description{
  This is the frontend for estimating Heckman-style selection models either with
  one or two outcomes (also known as generalized tobit models).
}

\usage{
selection(selection, outcome, data = sys.frame(sys.parent()),
   method = "ml", start = NULL, 
   ys = FALSE, xs = FALSE, yo = FALSE, xo = FALSE,
   mfs = FALSE, mfo = FALSE, print.level = 0, ...)

heckit( selection, outcome, data = sys.frame(sys.parent()),
   method = "2step", ... )
}

\arguments{
  \item{selection}{formula, the selection equation.}
  \item{outcome}{the outcome equation(s).  Either a single equation (for
  tobit 2 models), or a list of two equations (tobit 5 models).}
  \item{data}{an optional data frame, list or environment (or object
          coercible by \code{as.data.frame} to a data frame) containing the
          variables in the model.  If not found in \code{data}, the
          variables are taken from \code{environment(formula)}, typically
          the environment from which \code{selection} is called.}
  \item{method}{how to estimate the model.  Either "ml" for Maximum
    Likelihood, or "2step" for 2-step estimation.}
  \item{start}{initial values for the ML estimation.  Should have names.}
  \item{ys, yo, xs, xo, mfs, mfo}{logicals.  If true, the response (\code{y}),
    model matrix (\code{x}) or the model frame (\code{mf})
    of the selection (\code{s}) or outcome
    (\code{o}) equation(s) are returned.}
  \item{print.level}{integer.  Various debugging information, higher
    value gives more information.}
  \item{\dots}{additional parameters for the corresponding fitting
    functions \code{\link{tobit2fit}}, \code{\link{tobit5fit}},
    \code{\link{heckit2fit}}, and \code{\link{heckit5fit}}.}
}
\details{
   The endogenous variable of the argument 'selection' must have exactly
   two levels (e.g. 'FALSE' and 'TRUE', or '0' and '1').
   By default the levels are sorted in increasing order
   ('FALSE' is before 'TRUE', and '0' is before '1').

   For tobit-2 (sample selection) models,
   only those observations are included in the second step estimation
   (argument 'outcome'), where this variable equals the second element
   of its levels (e.g. 'TRUE' or '1').

   For tobit-5 (switching regression) models,
   in the second step the first outcome equation
   (first element of argument 'outcome') is estimated
   only for those observations,
   where this endogenous variable of the selections equation
   equals the first element of its levels (e.g. 'FALSE' or '0').
   The second outcome equation is estimated
   only for those observations,
   where this variable equals the second element of its levels
   (e.g. 'TRUE' or '1').

  \code{NA}-s are allowed in the data.  These are ignored if the
  corresponding outcome is unobserved, otherwise observations which
  contain \code{NA} (either in selection or outcome) are
  removed.
  
  These selection models assume a known (multivariate normal)
  distribution
  of error terms.  Because of this, the instruments (exclusion
  restrictions) are not necessary.  However, if no instruments are
  supplied, the results are based solely on the assumption on
  multivariate normality.  This may or may not be an appropriate
  assumption for a particular problem.

  The (generic) function 'coef' ('\code{\link{coef.selection}}')
  can be used to extract the estimated coefficients.
  The (generic) function 'vcov' ('\code{\link{vcov.selection}}')
  can be used to extract the estimated variance covariance matrix 
  of the coefficients.
  The (generic) function 'print' ('\code{print.selection}')
  can be used to print a few results.
  The (generic) function 'summary' ('\code{\link{summary.selection}}')
  can be used to obtain and print detailed results.
}
\value{
  'selection' returns an object of class "selection".
  If the model estimated by Maximum Likelihood (argument method = "ml"),
  this object is a list, which has all the components of 'maxLik', and in
  addition the elements 'twoStep', 'init', 'param', termS, termO,
  and if requested 'ys', 'xs', 'yo', 'xo', 'mfs', and 'mfo'.
  If a tobit-2 (sample selection) model is estimated by the
  two-step method (argument method = "2step"),
  the returned object is list with components
  'probit', 'coefficients', 'param', 'vcov', 'lm', 'sigma', 'rho',
  'invMillsRatio', and 'imrDelta'.
  If a tobit-2 (sample selection) model is estimated by the
  two-step method (argument method = "2step"),
  the returned object is list with components
  'coefficients', 'vcov', 'probit', 'lm1', 'lm2', 'rho1', 'rho2',
  'sigma1', 'sigma2', 'termsS', 'termsO', 'param',
  and if requested 'ys', 'xs', 'yo', 'xo', 'mfs', and 'mfo'.

  \item{probit}{object of class 'probit' that contains the results of the
    1st step (probit estimation) (only for two-step estimations).}
  \item{twoStep}{(only if initial values not given)
    results of the 2-step estimation, used for initial values}
  \item{init}{initial values for ML estimation}
  \item{termsS, termsO}{terms for the selection and outcome equation}
  \item{ys, xs, yo, xo, mfs, mfo}{response, matrix and frame of the
    selection- and outcome equations (as a list of two for the latter).
    NULL, if not requested.}
  \item{coefficients}{estimated coefficients, the complete model.
    coefficient for the Inverse Mills ratio is treated as a parameter
    (\eqn{=\varrho \sigma}{= rho * sigma}).}
  \item{vcov}{variance covariance matrix of the estimated coefficients.}
  \item{param}{a list with following components: \code{index}, a list of
    numeric vectors: where in the \code{coef} the component are
    located;
    \code{oIntercept}, a logical: whether the outcome equation includes
    intercept;
    \code{N0, N1}, integer, number of observations with unobserved and
    observed outcomes;
    \code{nObs}, integer, number of valid observations;
    \code{nParam}, integer, number of the parameters in the model (not
    all are independent);
    \code{df}, integer, degrees of freedom.  Note this is not equal to
    \code{nObs - nParam} because of the parameters are not independent
    in all the cases.
  }
  \item{lm, lm1, lm2}{objects of class 'lm' that contain the results
    of the 2nd step estimation(s) of the outcome equation(s).
    Note: the standard errors of this
    estimation are biased, because they do not account for the
    estimation of \eqn{\gamma} in the 1st step estimation
    (the correct standard errors are returned by \code{summary} and they
    are contained in \code{vcov} component).}
  \item{sigma, sigma1, sigma2}{the standard error(s)
    of the error terms of the outcome equation(s).}
  \item{rho, rho1, rho2}{the estimated correlation coefficient(s)
    between the error term of the selection equation
    and the outcome equation(s).}
  \item{invMillsRatio}{the inverse Mills Ratios calculated from the results of
    the 1st step probit estimation.}
  \item{imrDelta}{the \eqn{\delta}s calculated from the inverse Mills Ratios and
    the results of the 1st step probit estimation.}
}

\references{
  Cameron, A. C. and Trivedi, P. K. (2005)
  \emph{Microeconometrics: Methods and Applications},
  Cambridge University Press.

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Heckman, J. (1976)
  The common structure of statistical models of
  truncation, sample selection and limited dependent
  variables and a simple estimator for such models,
  \emph{Annals of Economic and Social Measurement},
  5(4), p. 475-492.

  Johnston, J. and J. DiNardo (1997)
  \emph{Econometric Methods, Fourth Edition}, McGraw-Hill.

  Lee, L., G. Maddala and R. Trost (1980)
  Asymetric covariance matrices of two-stage probit and two-stage tobit
  methods for simultaneous equations models with selectivity.
  \emph{Econometrica}, 48, p. 491-503.

  Wooldridge, J. M. (2003)
  \emph{Introductory Econometrics: A Modern Approach, 2e}, Thomson South-Western.
}

\note{
  The 2-step estimate of 'rho' may be outside of the
  \eqn{[-1,1]}{[-1,1]} interval.  In that case the standard errors of
  invMillsRatio may be meaningless.
  }

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet \email{otoomet@ut.ee}}

\seealso{\code{\link{lm}}, \code{\link{glm}}, \code{\link{binomial}}}

\examples{
## Greene( 2003 ): example 22.8, page 786
data( Mroz87 )
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
# Two-step estimation
summary( heckit( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 ) )
# ML estimation
summary( selection( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 ) )

## Wooldridge( 2003 ): example 17.5, page 590
data( Mroz87 )
# Two-step estimation
summary( heckit( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87,
   method = "2step" ) )

## Cameron and Trivedi (2005): Section 16.6, page 553ff
data( RandHIE )
subsample <- RandHIE$year == 2 & !is.na( RandHIE$educdec )
selectEq <- binexp ~ logc + idp + lpi + fmde + physlm + disea +
   hlthg + hlthf + hlthp + linc + lfam + educdec + xage + female +
   child + fchild + black
outcomeEq <- lnmeddol ~ logc + idp + lpi + fmde + physlm + disea +
   hlthg + hlthf + hlthp + linc + lfam + educdec + xage + female +
   child + fchild + black
# ML estimation
cameron <- selection( selectEq, outcomeEq, data = RandHIE[ subsample, ] )
summary( cameron )

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

## tobit-5 example
N <- 500
   library(mvtnorm)
   vc <- diag(3)
   vc[lower.tri(vc)] <- c(0.9, 0.5, 0.6)
   vc[upper.tri(vc)] <- vc[lower.tri(vc)]
   eps <- rmvnorm(N, rep(0, 3), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo1 <- runif(N)
   yo1 <- xo1 + eps[,2]
   xo2 <- runif(N)
   yo2 <- xo2 + eps[,3]
   a <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2))
   summary(a)

## tobit2 example
   vc <- diag(2)
   vc[2,1] <- vc[1,2] <- -0.7
   eps <- rmvnorm(N, rep(0, 2), vc)
   xs <- runif(N)
   ys <- xs + eps[,1] > 0
   xo <- runif(N)
   yo <- (xo + eps[,2])*(ys > 0)
   a <- selection(ys~xs, yo ~xo)
   summary(a)
}

\keyword{models}
\keyword{regression}
