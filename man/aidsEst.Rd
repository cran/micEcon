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
      maxiterIL = 50, tolIL = 1e-5, alpha0 = 0, TX = FALSE, ... )
}

\arguments{
   \item{pNames}{a vector of strings containing the names of the prices.}
   \item{wNames}{a vector of strings containing the names of the expenditure
      shares.}
   \item{xtName}{a string containing the variable name of total expenditure.}
   \item{data}{a data frame containing the data.}
   \item{ivNames}{a vector of strings containing the names of instrumental
      variables.}
   \item{qNames}{an optional vector of strings containing the names of
      the quantities (just to label elasticities).}
   \item{method}{the method to estimate the aids (see details).}
   \item{hom}{logical. Should the homogeneity condition be imposed?}
   \item{sym}{logical. Should the symmetry condition be imposed?}
   \item{elaFormula}{the elasticity formula (see \code{\link{aidsEla}}).}
   \item{pxBase}{The base to calculate the LA-AIDS price indices
      (see \code{\link{aidsPx}}).}
   \item{estMethod}{estimation method (e.g. 'SUR' or '3SLS',
      see \code{\link[systemfit]{systemfit}}).}
   \item{maxiterIL}{maximum number of iterations of the
      'Iterated Linear Least Squares Estimation'.}
   \item{tolIL}{tolerance level of the 'Iterated Linear Least Squares
      Estimation'.}
   \item{alpha0}{the intercept of the translog price index (\eqn{\alpha_0}).}
   \item{TX}{logical. Method to impose homogeneity and symmetry restrictions:
      either via R.restr (default) or via TX
      (see \code{\link[systemfit]{systemfit}}).}
   \item{...}{arguments passed to \code{\link[systemfit]{systemfit}}.}
}

\details{
   At the moment two basic estimation methods are available:
   The 'Linear Approximate AIDS' (LA) and the 'Iterative Linear Least Squares
   Estimator' (IL) proposed by Blundell and Robin (1999).\cr
   The LA-AIDS can be estimated with
   \begin{itemize}
   \item Stone price index ('LA:S'),
   \item Stone price index with lagged shares ('LA:SL'),
   \item loglinear analogue to the Paasche price index ('LA:P'),
   \item loglinear analogue of the Laspeyres price index ('LA:L'), and
   \item Tornqvist price index ('LA:T').
   \end{itemize}

   The 'Iterative Linear Least Squares Estimator' (IL) needs starting
   values for the (translog) price index. The price index used to
   calculate the initial price index can be specified in the same way
   as for the LA-AIDS (e.g. 'IL:L')
}

\value{
   a list of class \code{aidsEst} containing following objects:
   \item{coef}{a list containing the vectors/matrix of the estimated
      coefficients (alpha, beta, and gamma).}
   \item{ela}{a list containing the elasticities
      (see \code{\link{aidsEla}}).}
   \item{r2}{\eqn{R^2}-values of all share equations.}
   \item{r2q}{\eqn{R^2}-values of the estimated quantities.}
   \item{wFitted}{fitted expenditure shares.}
   \item{wResid}{residuals of the expenditure shares.}
   \item{qObs}{observed quantities / quantitiy indices.}
   \item{qFitted}{fitted quantities / quantitiy indices.}
   \item{qResid}{residuals of the estimated quantities.}
   \item{iter}{iterations of SUR/3SLS estimation(s).
      If the AIDS is estimated by the 'Iterated Linear Least Squares
      Estimator' (ILLE):
      a vector containing the SUR/3SLS iterations at each iteration.}
   \item{iterIL}{number of iterations of the 'Iterated Linear Least Squares
      Estimation'.}
   \item{method}{the method used to estimate the aids (see details).}
   \item{lnp}{log of the price index used for estimation.}
   \item{hom}{logical. Was the homogeneity condition imposed?}
   \item{sym}{logical. Was the symmetry condition imposed?}
   \item{estMethod}{estimation method (see \code{\link[systemfit]{systemfit}}).}
   \item{rcovformula}{formula used to calculate the estimated residual
      covariance matrix (see \code{\link[systemfit]{systemfit}}).}
   \item{pMeans}{means of the prices.}
   \item{wMeans}{means of the expenditure shares.}
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Blundell, R. and J.M. Robin (1999)
   Estimationin Large and Disaggregated Demand Systems:
   An Estimator for Conditionally Linear Systems.
   \emph{Journal of Applied Econometrics}, 14, p. 209-232.
}

\seealso{\code{\link{aidsEla}}, \code{\link{aidsCalc}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   # Using data published in Blanciforti, Green & King (1986)
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   ## Repeating the demand analysis of Blanciforti, Green & King (1986)
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, method = "LA:SL", elaFormula = "Ch",
      maxiter = 1, rcovformula = 1, tol = 1e-7 )
   print( estResult )

   ## Repeating the evaluation of different elasticity formulas of
   ## Green & Alston (1990)
   pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   # AIDS
   estResultA <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86[ -1, ], maxiter = 1,
      elaFormula = "AIDS", rcovformula=1, tol=1e-7,
      method = "IL:L" )
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
