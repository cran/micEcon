\name{aidsEla}
\alias{aidsEla}
\title{ Elasticities of the AIDS model }

\description{
   Calculates the demand elasticities of an AIDS model.
}

\usage{
   aidsEla( coef, W, P = NULL, formula = "AIDS",
      qNames = NULL, pNames = NULL )
}

\arguments{
   \item{coef}{a list containing the coefficients alpha, beta and gamma.}
   \item{W}{a vector of the shares at which the elasticities should be calculated.}
   \item{P}{a vector of the prices at which the elasticities should be calculated.}
   \item{formula}{the elsticity formulas to be used (see details).}
   \item{qNames}{an optional vector of strings containing the names of
      the quantities to label elasticities.}
   \item{pNames}{an optional vector of strings containing the names of
      the prices to label elasticities.}
}

\details{
   At the moment only the elasticity formulas of the orginal AIDS ('AIDS'),
   the formula of Chalfant (1987) ('Ch') and the formula of
   Eales and Unnevehr (1988) ('EU') are implemented.
}

\value{
  a list containing following elements:
  \item{formula}{the elasticity formula used to calculate these elasticities}
  \item{exp}{vector of expenditure elasticities}
  \item{marshall}{matrix of Marshallian (uncompensated) price elasticities}
  \item{hicks}{matrix of Hicksian (compensated) price elasticities}
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
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   estResult <- aidsEst( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, method = "LA:L" )
   wMeans <- colMeans( Blanciforti86[ , c( "wFood1", "wFood2",
      "wFood3", "wFood4" ) ] )
   aidsEla( estResult$coef, wMeans, formula = "Ch" )
}

\keyword{ models }
