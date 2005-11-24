\name{snqProfitFixEla}
\alias{snqProfitFixEla}
\title{Fixed Factor Elasticities of SNQ Profit function}

\description{
   Calculates the Fixed Factor Elasticities of a Symmetric Normalized
   Quadratic (SNQ) profit function.
}

\usage{ snqProfitFixEla( delta, gamma, quant, fix, weights )}

\arguments{
   \item{delta}{matrix of estimated \eqn{\delta} coefficients.}
   \item{gamma}{matrix of estimated \eqn{\gamma} coefficients.}
   \item{quant}{vector of netput quantities at which the elasticities
      should be calculated.}
   \item{fix}{vector of quantities of fixed factors at which the elasticities
      should be calculated.}
   \item{weights}{vector of weights of prices used for normalization.}
}

\note{
   A fixed factor elasticity is defined as
   \deqn{E_{ij} = \frac{ \displaystyle \frac{ \partial q_i }{ q_i } }
   {  \displaystyle \frac{ \partial z_j }{ z_j } } =
   \frac{ \partial q_i }{ \partial z_j } \cdot \frac{ z_j }{ q_i } }
   Thus, e.g. \eqn{E_{ij}=0.5} means that if the quantity of fixed factor j
   (\eqn{z_j}) increases by 1\%, the quantity of netput i (\eqn{q_i}) will
   increase by 0.5\%.
}

\seealso{\code{\link{snqProfitEst}} and \code{\link{snqProfitEla}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   # just a stupid simple example
   snqProfitFixEla( matrix(1:6/6,3,2 ), matrix(4:1/4,2 ), c(1,1,1), c(1,1),
      c(0.4,0.3,0.3) )

   # now with real data
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput   / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time      <- c( 0:19 )
   pNames <- c( "pOutput", "pVarInput", "pLabor" )
   qNames <- c( "qOutput", "qVarInput", "qLabor" )
   fNames <- c( "land", "time" )

   estResult <- snqProfitEst( pNames, qNames, fNames, data=germanFarms )

   estResult$fixEla  # price elasticities at mean quantities of netputs
                     # and fixed factors

   # fixed factor elasticities at the last observation (1994/95)
   snqProfitFixEla( estResult$coef$delta, estResult$coef$gamma,
      estResult$estData[ 20, qNames ], estResult$estData[ 20, fNames ],
      estResult$weights )
}

\keyword{models}
