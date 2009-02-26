\name{snqProfitHessianDeriv}
\alias{snqProfitHessianDeriv}
\title{SNQ Profit function: Derivatives of the Hessian}

\description{
   Returns the matrix of derivatives of the vector of linear independent
   values of the Hessian with respect to the vector of the linear independent
   coefficients.
}

\usage{ snqProfitHessianDeriv( prices, weights, nFix = 0, form = 0 )}

\arguments{
   \item{prices}{vector of netput prices at which the derivatives
      should be calculated.}
   \item{weights}{vector of weights for normalizing prices.}
   \item{nFix}{number of (quasi-)fix inputs.}
   \item{form}{the functional form to be estimated (see
      \code{\link{snqProfitEst}}).}
}

\seealso{\code{\link{snqProfitHessian}}.}

\author{Arne Henningsen}

\examples{
   # just a stupid simple example
   snqProfitHessianDeriv( c(1,2,3),c(0.4,0.3,0.3) )

   # now with real data
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput   / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time <- c( 0:19 )
   priceNames <- c( "pOutput", "pVarInput", "pLabor" )
   quantNames <- c( "qOutput", "qVarInput", "qLabor" )

   estResult <- snqProfitEst( priceNames, quantNames, c("land","time"), data=germanFarms )

   snqProfitHessianDeriv( estResult$pMean, estResult$weights, 2 )
}

\keyword{models}
