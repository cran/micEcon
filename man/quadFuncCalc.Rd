\name{quadFuncCalc}
\alias{quadFuncCalc}
\title{Calculate dependent variable of a quadratic function}

\description{
   Calculate the dependent variable of a quadratic function.
}

\usage{
   quadFuncCalc( xNames, data, allCoef, quadHalf = TRUE  )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{allCoef}{vector containing all coefficients.}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
}

\value{
   a vector containing the endogenous variable.
}

\seealso{\code{\link{quadFuncEst}} and \code{\link{quadFuncDeriv}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # estimate a quadratic production function
   estResult <- quadFuncEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms )

   quadFuncCalc( c( "qLabor", "land", "qVarInput", "time" ), germanFarms,
      estResult$allCoef )
   #equal to estResult$fitted
}

\keyword{models}
