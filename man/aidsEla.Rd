\name{aidsEla}
\alias{aidsEla}
\title{ Elasticities of the AIDS model }

\description{
   Calculates the demand elasticities of an AIDS model.
}

\usage{
   aidsEla( coef, W, P = NULL, formula = "AIDS",
      qNames = NULL, pNames = NULL, coefVcov = NULL, df = NULL )
}

\arguments{
   \item{coef}{a list containing the coefficients alpha, beta and gamma.}
   \item{W}{a vector of the shares at which the elasticities should be calculated.}
   \item{P}{a vector of the prices at which the elasticities should be calculated
      (only needed if formula 'AIDS' is used).}
   \item{formula}{the elsticity formula to be used (see details).}
   \item{qNames}{an optional vector of strings containing the names of
      the quantities to label elasticities.}
   \item{pNames}{an optional vector of strings containing the names of
      the prices to label elasticities.}
   \item{coefVcov}{variance covariance matrix of the coefficients (optional).}
   \item{df}{degrees of freedom to calculate P-values of the elasticities 
      (optional).}
}

\details{
   At the moment only the elasticity formulas of the orginal AIDS ('AIDS'),
   the formula of Chalfant (1987) ('Ch') and the formula of
   Eales and Unnevehr (1988) ('EU') are implemented.
   The variance covariance matrices of the elasticities are calculated using
   the formula of Klein (1953, p. 258) (also known as the delta method).
   At the moment this is implemented only for the elasticity formulas of the 
   orginal AIDS.
}

\value{
  a list of class \code{aidsEla} containing following elements:
  \item{formula}{the elasticity formula used to calculate these elasticities.}
  \item{exp}{vector of expenditure elasticities.}
  \item{hicks}{matrix of Hicksian (compensated) price elasticities.}
  \item{marshall}{matrix of Marshallian (uncompensated) price elasticities.}
  \item{allVcov}{variance covariance matrix of all elasticities.}
  \item{expVcov}{variance covariance matrix of the expenditure elasticities.}
  \item{hicksVcov}{variance covariance matrix of the Hicksian (compensated) 
    price elasticities.}
  \item{marshallVcov}{variance covariance matrix of the Marshallian 
    (uncompensated) price elasticities.}
  \item{expStEr}{standard errors of the expenditure elasticities.}
  \item{hicksStEr}{standard errors of the Hicksian (compensated) price 
    elasticities.}
  \item{marshallStEr}{standard errors of the Marshallian (uncompensated) 
    price elasticities.}
  \item{expTval}{t-values of the expenditure elasticities.}
  \item{hicksTval}{t-values of the Hicksian (compensated) price elasticities.}
  \item{marshallTval}{t-values of the Marshallian (uncompensated) price 
    elasticities.}
  \item{expPval}{P-values of the expenditure elasticities.}
  \item{hicksPval}{P-values of the Hicksian (compensated) price elasticities.}
  \item{marshallPval}{P-values of the Marshallian (uncompensated) price 
    elasticities.}
}

\references{
   Chalfant, J.A. (1987)
   A Globally Flexible, Almost Ideal Demand System.
   \emph{Journal of Business and Economic Statistics}, 5, p. 233-242.

   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Eales J.S. and L.J. Unnevehr (1988)
   Demand for beef and chicken products: separability and structural change.
   \emph{American Journal of Agricultural Economics}, 70, p. 521-532.
   
   Klein L.R. (1953)
   \emph{A Textbook of Econometrics}. Row, Petersen and Co., New York.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, method = "LA:L" )
   wMeans <- colMeans( Blanciforti86[ , c( "wFood1", "wFood2",
      "wFood3", "wFood4" ) ] )
   aidsEla( estResult$coef, wMeans, formula = "Ch" )
}

\keyword{ models }
