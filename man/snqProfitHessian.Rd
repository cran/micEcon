\name{snqProfitHessian}
\alias{snqProfitHessian}
\title{SNQ Profit function: Hessian matrix}
\description{
   Returns the Hessian (substitution) matrix of a Symmetric Normalized
   Quadratic (SNQ) Profit Function.
}

\usage{snqProfitHessian( beta, prices, weights,
   scalingFactors = rep( 1, length( weights ) ) )}

\arguments{
   \item{beta}{matrix of the \eqn{beta} coefficients.}
   \item{prices}{vector of netput prices at which the Hessian
      should be calculated.}
   \item{weights}{vector of weights of prices for normalization.}
   \item{scalingFactors}{factors to scale prices (and quantities).}
}

\seealso{\code{\link{snqProfitEst}}, \code{\link{snqProfitEla}} and
   \code{\link{snqProfitHessianDeriv}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   # just a stupid simple example
   snqProfitHessian( matrix(101:109,3,3), c(1,1,1), c(0.4,0.3,0.3) )

   # now with real data
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput   / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time <- c( 0:19 )
   pNames <- c( "pOutput", "pVarInput", "pLabor" )
   qNames <- c( "qOutput", "qVarInput", "qLabor" )

   estResult <- snqProfitEst( pNames, qNames, c("land","time"), data=germanFarms )

   estResult$hessian  # the Hessian at mean prices and mean quantities

   # Hessian at the last observation (1994/95)
   snqProfitHessian( estResult$coef$beta, estResult$data[ 20, pNames ],
      estResult$weights, estResult$scalingFactors )
}

\keyword{models}
