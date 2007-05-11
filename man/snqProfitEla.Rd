\name{snqProfitEla}
\alias{snqProfitEla}
\title{Price Elasticities of SNQ Profit function}

\description{
   Calculates the Price Elasticities of a Symmetric Normalized Quadratic (SNQ)
   profit function.
}

\usage{ snqProfitEla(  beta, prices, quant, weights,
   scalingFactors = rep( 1, length( weights ) ),
   coefVcov = NULL, df = NULL )}

\arguments{
   \item{beta}{matrix of estimated \eqn{\beta} coefficients.}
   \item{prices}{vector of netput prices at which the elasticities
      should be calculated.}
   \item{quant}{vector of netput quantities at which the elasticities
      should be calculated.}
   \item{weights}{vector of weights of prices used for normalization.}
   \item{scalingFactors}{factors to scale prices (and quantities).}
   \item{coefVcov}{variance covariance matrix of the coefficients (optional).}
   \item{df}{degrees of freedom to calculate P-values of the elasticities
      (optional).}
}

\note{
   A price elasticity is defined as
   \deqn{E_{ij} = \frac{ \displaystyle \frac{ \partial q_i }{ q_i } }
   {  \displaystyle \frac{ \partial p_j }{ p_j } } =
   \frac{ \partial q_i }{ \partial p_j } \cdot \frac{ p_j }{ q_i } }
   Thus, e.g. \eqn{E_{ij}=0.5} means that if the price of netput j (\eqn{p_j})
   increases by 1\%, the quantity of netput i (\eqn{q_i}) will
   increase by 0.5\%.
}

\value{
  a list of class \code{snqProfitEla} containing following elements:
  \item{ela}{matrix of the price elasticities.}
  \item{vcov}{variance covariance matrix of the price elasticities.}
  \item{stEr}{standard errors of the price elasticities.}
  \item{tval}{t-values of the price elasticities.}
  \item{pval}{P-values of the price elasticities.}
}

\seealso{\code{\link{snqProfitEst}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   # just a stupid simple example
   snqProfitEla( matrix(101:109,3,3), c(1,1,1), c(1,-1,-1), c(0.4,0.3,0.3) )

   # now with real data
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput   / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time      <- c( 0:19 )
   priceNames <- c( "pOutput", "pVarInput", "pLabor" )
   quantNames <- c( "qOutput", "qVarInput", "qLabor" )

   estResult <- snqProfitEst( priceNames, quantNames, c("land","time"), data=germanFarms )

   estResult$ela  # price elasticities at mean prices and mean quantities

   # price elasticities at the last observation (1994/95)
   snqProfitEla( estResult$coef$beta, estResult$data[ 20, priceNames ],
      estResult$data[ 20, quantNames ], estResult$weights,
      estResult$scalingFactors )
}

\keyword{models}
