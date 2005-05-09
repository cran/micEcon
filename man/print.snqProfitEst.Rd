\name{print.snqProfitEst}
\alias{print.snqProfitEst}
\title{Print output of estimated SNQ profit function}

\description{
  This function prints a summary estimation results of
  a symmetric normalized quadraric (SNQ) profit function.
}

\usage{
   \method{print}{snqProfitEst}( x, \dots )
}

\arguments{
  \item{x}{an object of class \code{snqProfitEst}.}
  \item{\dots}{arguments passed to \code{\link{print}}.}
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{snqProfitEst}}}

\examples{
\dontrun{library( systemfit )}
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   pNames <- c( "pOutput", "pVarInput", "pLabor" )
   qNames <- c( "qOutput", "qVarInput", "qLabor" )

   estResult <- snqProfitEst( pNames, qNames, "land", data = germanFarms )
   print( estResult )
}

\keyword{models}



