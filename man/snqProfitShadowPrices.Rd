\name{snqProfitShadowPrices}
\alias{snqProfitShadowPrices}
\title{Shadow Prices of a SNQ Profit function}

\description{
   Calculates the shadow prices of a Symmetric Normalized Quadratic (SNQ)
   profit function.
}

\usage{
   snqProfitShadowPrices( pNames, fNames, estResult = NULL,
   data = estResult$estData, weights = estResult$weights,
   coef = estResult$coef, form = estResult$form )
}

\arguments{
   \item{pNames}{a vector of strings containing the names of netput prices.}
   \item{fNames}{an optional vector of strings containing the names of the
      quantities of (quasi-)fix inputs.}
   \item{estResult}{object returned by \code{\link{snqProfitEst}}.}
   \item{data}{a data frame containing the data.}
   \item{weights}{vector of weights of prices used for normalization.}
   \item{coef}{a list containing the coefficients (at least delta and gamma).}
   \item{form}{the functional form to be estimated (see details).}
}

\seealso{
   \code{\link{snqProfitEst}}, \code{\link{snqProfitCalc}} and
   \code{\link{snqProfitEla}}.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput   / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time      <- c( 0:19 )
   pNames <- c( "pOutput", "pVarInput", "pLabor" )
   qNames <- c( "qOutput", "qVarInput", "qLabor" )
   fNames <- c( "land", "time" )

   estResult <- snqProfitEst( pNames, qNames, fNames, data = germanFarms )

   snqProfitShadowPrices( pNames, fNames, estResult )
}

\keyword{models}
