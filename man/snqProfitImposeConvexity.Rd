\name{snqProfitImposeConvexity}
\alias{snqProfitImposeConvexity}
\title{Imposing Convexity on a SNQ Profit function}
\description{
   Imposing Convexity on a Symmetric Normalized Quadratic (SNQ) Profit function.
}

\usage{snqProfitImposeConvexity( estResult, rankReduction = 0,
   start = 10, ... )
}

\arguments{
   \item{estResult}{object returned by \code{\link{snqProfitEst}}.}
   \item{rankReduction}{an integer specifying the reduction of the rank
      of the \eqn{\beta} matrix.}
   \item{start}{starting values of the triangular Cholesky matrix.}
   \item{...}{arguments passed to \code{\link{optim}}}
}

\details{
   The procedure proposed by Koebel, Falk and Laisney (2000, 2003) is
   applied to impose convexity in prices on an estimated symmetric
   normalized quadratic (SNQ) profit function.
}

\value{
   a list of class \code{snqProfitImposeConvexity} containing following objects:
   \item{mindist}{object returned by \code{\optim}.}
   \item{coef}{a list containing the vectors/matrix of the estimated
      coefficients:\cr
      * alpha = \eqn{\alpha_i}.\cr
      * beta = \eqn{\beta_{ij}}.\cr
      * delta =  \eqn{\delta_{ij}} (only if quasi-fix inputs are present).\cr
      * gamma = \eqn{\gamma_{ij}} (only if quasi-fix inputs are present).\cr
      * allCoef = vector of all coefficients.\cr
   }
   \item{ela}{matrix with the price elasticities at mean prices and mean quantities.}
   \item{hessian}{hessian matrix of the profit function with respect to prices
      evaluated at mean prices.}
   \item{convexity}{logical. Convexity of the profit function.}
   \item{r2}{\eqn{R^2}-values of all netput equations.}
   \item{weights}{the weights of prices used for normalization.}
   \item{normPrice}{vector used for normalization of prices.}
   \item{estData}{data frame used for estimation (contains the
      (scaled) netput prices, (scaled) netput quantities,
      (not scaled) fix inputs and the price index used for normalization.}
   \item{fitted}{data frame that contains the fitted netput quantities and
      the fitted profit.}
   \item{residuals}{data frame that contains the residuals of the netput
      quantities.}
   \item{form}{the functional form (see \code{\link{snqProfitEst}}).}
   \item{pMeans}{means of the (scaled) netput prices.}
   \item{qMeans}{means of the (scaled) netput quantities.}
   \item{fMeans}{means of the (quasi-)fix input quantities.}
}

\references{
   Koebel, B., M. Falk and F. Laisney (2000),
   Imposing and Testing Curvature Conditions on a Box-Cox Cost Function.
   Discussion Paper No. 00-70, ZEW, Mannheim,
   \url{ftp://ftp.zew.de/pub/zew-docs/dp/dp0070.pdf}.

   Koebel, B., M. Falk and F. Laisney (2003),
   Imposing and Testing Curvature Conditions on a Box-Cox Cost Function.
   \emph{Journal of Business and Economic Statistics}, 21, p. 319-335.
}

\seealso{\code{\link{snqProfitEst}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time <- c( 0:19 )
   pNames <- c( "pOutput", "pVarInput", "pLabor" )
   qNames <- c( "qOutput", "qVarInput", "qLabor" )
   fNames <- c( "land", "time" )
   estResult <- snqProfitEst( pNames, qNames, fNames, data = germanFarms )
   estResult  # Note: it is NOT convex in netput prices
   estResultConvex <- snqProfitImposeConvexity( estResult )
   estResultConvex  # now it is convex
}

\keyword{models}

