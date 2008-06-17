\name{translogCalc}
\alias{translogCalc}
\title{Calculate dependent variable of a translog function}

\description{
   Calculate the dependent variable of a translog function.
}

\usage{
   translogCalc( xNames, data, coef, quadHalf = TRUE, dataLogged = FALSE  )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients.}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   A vector containing the endogenous variable.
   If the inputs are provided as logarithmic values
   (argument \code{dataLogged} is \code{TRUE}),
   the endogenous variable is returned as logarithm;
   non-logarithmic values are returned otherwise.
}

\seealso{\code{\link{translogEst}} and \code{\link{translogDeriv}}.}

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
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms )

   translogCalc( c( "qLabor", "land", "qVarInput", "time" ), germanFarms,
      coef( estResult ) )
   #equal to estResult$fitted
}

\keyword{models}
