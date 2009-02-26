\name{translogCalc}
\alias{translogCalc}
\title{Calculate dependent variable of a translog function}

\description{
   Calculate the dependent variable of a translog function.
}

\usage{
translogCalc( xNames, data, coef, shifterNames = NULL,
   quadHalf = TRUE, dataLogged = FALSE  )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients: 
      if there are \code{n} exogenous variables in \code{xNames}
      and \code{m} shifter variables in \code{shifterNames},
      the \code{n+1} alpha coefficients must have names 
      \code{a_0}, \ldots, \code{a_n},
      the \code{n*(n+1)/2} beta coefficients must have names
      \code{b_1_1}, \ldots, \code{b_1_n}, \ldots, \code{b_n_n},
      and the \code{m} delta coefficients must have names
      \code{d_1}, \ldots, \code{d_m}
      (only the elements of the upper right triangle of the beta matrix
      are directly obtained from \code{coef};
      the elements of the lower left triangle are obtained by assuming
      symmetry of the beta matrix).}
   \item{shifterNames}{a vector of strings containing the names
      of the independent variables that should be included
      as shifters only (not in quadratic or interaction terms).}
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

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # estimate a Translog production function
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms )

   translogCalc( c( "qLabor", "land", "qVarInput", "time" ), germanFarms,
      coef( estResult ) )
   #equal to estResult$fitted
}

\keyword{models}
