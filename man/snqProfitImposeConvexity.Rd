\name{snqProfitImposeConvexity}
\alias{snqProfitImposeConvexity}
\title{Imposing Convexity on a SNQ Profit function}
\description{
   Imposing Convexity on a Symmetric Normalized Quadratic (SNQ) Profit function.
}

\usage{snqProfitImposeConvexity( estResult, rankReduction = 0,
   start = 10, optimMethod = "BFGS", control = list( maxit=5000 ),
   stErMethod = "none", nRep = 1000, verbose = 0 )
}

\arguments{
   \item{estResult}{object returned by \code{\link{snqProfitEst}}.}
   \item{rankReduction}{an integer specifying the reduction of the rank
      of the \eqn{\beta} matrix.}
   \item{start}{starting values of the triangular Cholesky matrix.}
   \item{optimMethod}{method to be used by \code{\link{optim}}.}
   \item{control}{list of control parameters passed to \code{\link{optim}}.}
   \item{stErMethod}{method to compute standard errors, either
      'none', 'resample', 'jackknife' or 'coefSim' (see details).}
   \item{nRep}{number of replications to compute the standard errors
      if \code{stErMethod} is either 'resample' or 'coefSim'.}
   \item{verbose}{an integer idicating the verbose level.}
}

\details{
   The procedure proposed by Koebel, Falk and Laisney (2000, 2003) is
   applied to impose convexity in prices on an estimated symmetric
   normalized quadratic (SNQ) profit function.\cr
   The standard errors of the restricted coefficients can be either
   calculated by
   bootstrap resampling ('resampling'),
   jackknife ('jacknife') or
   by simulating the distribution of the unrestricted coefficients
   using its variance covariance matrix ('coefSim').
}

\value{
   a list of class \code{snqProfitImposeConvexity} containing the same
   objects as an object of class \code{\link{snqProfitEst}}
   and additionally the objects:
   \item{mindist}{object returned by \code{\link{optim}}.}
   \item{sim}{results of the simulation to obtain the standard errors
      of the estimated coefficients.}
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

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   priceNames <- c( "pOutput", "pVarInput", "pLabor" )
   quantNames <- c( "qOutput", "qVarInput", "qLabor" )
   estResult <- snqProfitEst( priceNames, quantNames, "land", data = germanFarms )
   estResult  # Note: it is NOT convex in netput prices
   estResultConvex <- snqProfitImposeConvexity( estResult )
   estResultConvex  # now it is convex
}

\keyword{models}

