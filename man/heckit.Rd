\name{heckit}
\alias{heckit}

\title{2-step Heckman (heckit) estimation}

\description{
   \code{heckit} performs a 2-step Heckman (heckit) estimation
   that corrects for non-random sample selection.
}

\usage{
   heckit( formula, probitformula, data, inst = NULL, print.level = 0 )
}

\arguments{
  \item{formula}{formula to be estimated}
  \item{probitformula}{formula for the probit estimation (1st step)}
  \item{data}{a data frame containing the variables in the model}
  \item{inst}{one-sided formula specifying instrumental variables for a
    2SLS/IV estimation on the second step.}
  \item{print.level}{this argument determines the level of printing which is
    done during the minimization process.  The default value of
    '0' means that no printing occurs, a value of '1' means that
    heckit reports what is currently done.}
}

\value{
   \code{heckit} returns an object of class 'heckit' containing
   following elements:

   \item{coef}{estimated coefficients, standard errors, t-values and
   p-values}

   \item{vcov}{variance covariance matrix of the estimated coefficients}

   \item{probit}{object of class 'glm' that contains the results of the
   1st step (probit estimation).}

   \item{lm}{object of class 'lm' that contains the results of the
   2nd step (linear estimation). Note: the standard errors of this
   estimation are biased, because they do not account for the
   estimation of \eqn{\gamma} in the 1st step estimation
   (the correct standard errors are returned in \code{coef}}.

   \item{sigma}{the estimated \eqn{\sigma}, the standard error of the residuals.}

   \item{rho}{the estimated \eqn{\rho}, see Greene (2003, p. 784).}

   \item{probitLambda}{the \eqn{\lambda}s based on the results of the
   1sr step probit estimation (also known as inverse Mills ratio).}

   \item{probitDelta}{the \eqn{\delta}s based on the results of the
   1sr step probit estimation.}
}

\references{

  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Johnston, J. and J. DiNardo (1997)
  \emph{Econometric Methods, Fourth Edition}, McGraw-Hill.

  Wooldridge, J. M. (2003)
  \emph{Introductory Econometrics: A Modern Approach, 2e}, Thomson South-Western.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{lm}}, \code{\link{glm}}, \code{\link{binomial}}}


\examples{
## Greene( 2003 ): example 22.8, page 786
data( Mroz87 )
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
greene <- heckit( wage ~ exper + I( exper^2 ) + educ + city,
   lfp ~ age + I( age^2 ) + faminc + kids + educ, Mroz87 )
summary( greene )        # print summary
summary( greene$probit ) # summary of the 1st step probit estimation
                         # this is Example 21.4, p. 681f
greene$sigma             # estimated sigma
greene$rho               # estimated rho

## Wooldridge( 2003 ): example 17.5, page 590
data( Mroz87 )
wooldridge <- heckit( log( wage ) ~ educ + exper + I( exper^2 ),
   lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age + kids5 + kids618, Mroz87 )
summary( wooldridge )        # summary of the 1st step probit estimation
                             # (Example 17.1, p. 562f) and 2nd step OLS regression
wooldridge$sigma             # estimated sigma
wooldridge$rho               # estimated rho

## example using random numbers
nObs <- 1000
myData <- data.frame( no = c( 1:nObs ), x1 = rnorm( nObs ), x2 = rnorm( nObs ) )
myData$y <- 2 + myData$x1 + 0.9 * rnorm( nObs )
myData$s <- ( 2 * myData$x1 + myData$x2 + 4 * rnorm( nObs ) - 0.2 ) > 0
myData$y[ !myData$s ] <- NA
myHeckit <- heckit( y ~ x1, s ~ x1 + x2, myData, print.level = 1 )

## example using random numbers with IV/2SLS estimation
nObs <- 1000
myData <- data.frame( no = c( 1:nObs ), x1 = rnorm( nObs ), x2 = rnorm( nObs ),
   u = 0.5 * rnorm( nObs ) )
myData$w <- 1 + myData$x1 + 0.2 * myData$u + 0.1 * rnorm( nObs )
myData$y <- 2 + myData$w + myData$u
myData$s <- ( 2 * myData$x1 + myData$x2 + 4 * rnorm( nObs ) - 0.2 ) > 0
myData$y[ !myData$s ] <- NA
myHeckit <- heckit( y ~ w, s ~ x1 + x2, data = myData,
   inst = ~ x1, print.level = 1 )
}

\keyword{models}
