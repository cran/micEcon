\name{aidsEst}
\alias{aidsEst}
\title{Estimation of the Almost Ideal Demand System (AIDS) }

\description{
   \code{aidsEst} does a full demand analysis with the Almost Ideal
   Demand System (AIDS): econometric estimation, calculation
   of elasticities, ...
}

\usage{ aidsEst( pNames, wNames, xtName, data = NULL,
      ivNames = NULL, qNames = wNames,
      method = "LA:L", hom = TRUE, sym = TRUE,
      elaFormula = "Ch", pxBase,
      estMethod = ifelse( is.null( ivNames ), "SUR", "3SLS" ),
      maxiterMk = 50, tolMk = 1e-5, alpha0 = 0, ... )
}

\arguments{
   \item{pNames}{a vector of strings containing the names of prices.}
   \item{wNames}{a vector of strings containing the names of expenditure
      shares.}
   \item{xtName}{a name of the total expenditure variable.}
   \item{data}{a data frame containing the data.}
   \item{ivNames}{a vector of strings containing the names of instrumental
      variables.}
   \item{qNames}{an optional vector of strings containing the names of
      the quantities (just to label elasticities).}
   \item{method}{the method to estimate the aids (see details).}
   \item{hom}{should the homogeneity condition be imposed?}
   \item{sym}{should the symmetry condition be imposed?}
   \item{elaFormula}{the elasticity formula (see \code{\link{aidsEla}}).}
   \item{pxBase}{The base to calculate the LA-AIDS price indices
      (see \code{\link{aidsPx}}).}
   \item{estMethod}{estimation method (e.g. "SUR",
      see \code{\link[systemfit]{systemfit}}).}
   \item{maxiterMk}{maximum number of iterations of the
      Michalek and Keyzer estimation}
   \item{tolMk}{tolerance level of the Michalek and Keyzer iteration}
   \item{alpha0}{the intercept of the translog price index (\eqn{\alpha_0})}
   \item{...}{argunents passed to \code{\link[systemfit]{systemfit}}}
}

\details{
   At the moment two basic estimation methods are available:
   The 'Linear Approximate AIDS' (LA) and the iterative estimation method
   of Michalek and Keyzer (MK).\cr
   The LA-AIDS can be estimated with\cr
   * Stone price index ('LA:S'),\cr
   * Stone price index with lagged shares ('LA:SL'),\cr
   * loglinear analogue to the Paasche price index ('LA:P'),\cr
   * loglinear analogue of the Laspeyres price index ('LA:L'), and\cr
   * Tornqvist price index ('LA:T').

   The iterative estimation method of Michalek and Keyzer needs starting
   values for the (translog) price index. Which price index is used to
   calculate the initial price index can be specified in the same way
   as for the LA-AIDS (e.g. 'MK:L')
}

\value{
   a list of class \code{aidsEst} containing following objects:
   \item{coef}{a list containing the vectors/matrix of the estimated
      coefficients (alpha = \eqn{\alpha_i}, beta = \eqn{\beta_i},
      gamma = \eqn{\gamma_{ij}}.}
   \item{ela}{a list containing the elasticities
      (see \code{\link{aidsEla}}).}
   \item{r2}{\eqn{R^2}-values of all budget share equations.}
   \item{r2q}{\eqn{R^2}-values of the estimated quantities.}
   \item{wFitted}{fitted budget shares.}
   \item{wResid}{residuals of the budget shares.}
   \item{qObs}{observed quantities / quantitiy indices.}
   \item{qFitted}{fitted quantities / quantitiy indices.}
   \item{qResid}{residuals of the estimated quantities.}
   \item{iter}{iterations of SUR/3SLS estimation a vector containing
      the iterations of each iteration of the Michalek and Keyzer
      estimation.}
   \item{iterMk}{iterations of the Michalek and Keyzer estimation.}
   \item{method}{the method used to estimate the aids (see details).}
   \item{lnp}{log of the price index used for estimation.}
   \item{hom}{was the homogeneity condition imposed?}
   \item{sym}{was the symmetry condition imposed?}
   \item{estMethod}{estimation method (see \code{\link[systemfit]{systemfit}}).}
   \item{rcovformula}{formula used to calculate the estimated residual
      covariance matrix (see \code{\link[systemfit]{systemfit}}).}
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Michalek, J. and M. A. Keyzer (1992)
   Estimation of a two-stage LES-AIDS consumer demand system
   for eight EC countries.
   \emph{European Review of Agricultural Economics}, 19 (2), p. 137-163.
}

\seealso{\code{\link{aidsEla}}, \code{\link{aidsCalc}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   ## Repeating the demand analysis of Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, method = "LA:SL", elaFormula = "Ch",
      maxiter = 1, rcovformula = 1, tol = 1e-7 )
   print( estResult )

   ## Repeating the evaluation of different elasticity formulas of
   ## Green & Alston (1990)
   data( Blanciforti86 )
   pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   # AIDS
   estResultA <- aidsEst( pNames, wNames, "xFood",
    data = Blanciforti86[ 2:nrow( Blanciforti86 ), ], maxiter = 1,
      elaFormula = "AIDS", rcovformula=1, tol=1e-7,
   method = "MK:L" )
   print( diag( estResultA$ela$marshall ) )

   # LA-AIDS + formula of AIDS
   estResultL1 <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86, maxiter = 1, elaFormula = "AIDS",
      rcovformula=1, tol=1e-7, method = "LA:SL" )
   print( diag( estResultL1$ela$marshall ) )

   # LA-AIDS + formula of Eales + Unnevehr
   estResultL2 <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86, maxiter = 1, elaFormula = "EU",
      rcovformula=1, tol=1e-7, method = "LA:SL" )
   print( diag( estResultL2$ela$marshall ) )

   # LA-AIDS + formula of Chalfant:
   estResultL3 <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86, maxiter = 1, elaFormula = "Ch",
      rcovformula=1, tol=1e-7, method = "LA:SL" )
   print( diag( estResultL3$ela$marshall ) )
}

\keyword{models}
