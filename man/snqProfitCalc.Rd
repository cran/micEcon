\name{snqProfitCalc}
\alias{snqProfitCalc}
\title{Calculations with the SNQ Profit function}

\description{
   Calculation of netput quantities and profit with the
   Symmetric Normalized Quadratic (SNQ) Profit function.
}

\usage{
   snqProfitCalc( pNames, fNames, data, weights, coef, form = 0  )
}

\arguments{
   \item{pNames}{a vector of strings containing the names of netput prices.}
   \item{fNames}{an optional vector of strings containing the names of the
      quantities of (quasi-)fix inputs.}
   \item{data}{a data frame containing the data.}
   \item{weights}{vector of weights of the prices for normalization.}
   \item{coef}{a list containing the coefficients alpha, beta, delta and gamma.}
   \item{form}{the functional form to be estimated (see \code{\link{snqProfitEst}}).}
}


\value{
   a data frame: the first n columns are the netput quantities,
   the last column is the profit.
}

\references{
   Diewert, W.E. and T.J. Wales (1987)
   Flexible functional forms and global curvature conditions.
   \emph{Econometrica}, 55, p. 43-68.

   Diewert, W.E. and T.J. Wales (1992)
   Quadratic Spline Models for Producer's Supply and Demand Functions.
   \emph{International Economic Review}, 33, p. 705-722.

   Kohli, U.R. (1993)
   A symmetric normalized quadratic GNP function and the US demand
   for imports and supply of exports.
   \emph{International Economic Review}, 34, p. 243-255.
}

\seealso{\code{\link{snqProfitEst}} and \code{\link{snqProfitWeights}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}


\examples{
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time      <- c( 0:19 )
   pNames <- c( "pOutput", "pVarInput", "pLabor" )
   qNames <- c( "qOutput", "qVarInput", "qLabor" )
   fNames <- c( "land", "time" )

   estResult <- snqProfitEst( pNames, qNames, fNames, data = germanFarms )
   snqProfitCalc( pNames, fNames, estResult$estData, estResult$weights,
      estResult$coef )
}

\keyword{models}

