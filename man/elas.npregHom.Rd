\name{elas.npregHom}
\alias{elas.npregHom}
\title{Elasticities of a Homogeneous Nonparametric Function}

\description{
   Calculate elasticities from a nonparametric regression
   with homogeneity of degree zero imposed on some variables.
}

\usage{
\method{elas}{npregHom}( object, data = NULL, yObs = FALSE, ... )
}

\arguments{
   \item{object}{object of class \code{npregHom}
      (returned by \code{\link{npregHom}}).}
   \item{data}{dataframe or a vector with named elements containing the data;
      if it is not specified,
      the data frame that was used for the nonparametric estimation
      is used for calculating elasticities.}
   \item{yObs}{logical. Use observed values of the andogenous variable.
      If \code{FALSE} (default) fitted values are used.}
   \item{\ldots}{currently ignored.}
}

\value{
   A data.frame,
   where each column corresponds to one of the continuous independent variables.
}

\seealso{
   \code{\link{npregHom}}.
}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # weights to impose
   weights <- c(
      pOutput = mean( germanFarms$qOutput ),
      pVarInput = mean( germanFarms$qVarInput ),
      pLabor = mean( germanFarms$qLabor ) )
   weights <- weights / sum( weights )

   # estimate an input demand function
   estResult <- npregHom( "qVarInput",
      xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
      data = germanFarms, homWeights = weights )

   # calculate elasticities
   elas( estResult )
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
