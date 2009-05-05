\name{cobbDouglasCalc}
\alias{cobbDouglasCalc}
\title{Calculate dependent variable of a Cobb-Douglas function}

\description{
   Calculate the dependent variable of a Cobb-Douglas function.
}

\usage{
cobbDouglasCalc( xNames, data, coef, coefCov = NULL, dataLogged = FALSE  )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{data frame containing the data.}
   \item{coef}{vector containing the coefficients:
      if the elements of the vector have no names,
      the first element is taken as intercept of the \emph{logged} equation
      and the following elements are taken as coefficients of
      the independent variables defined in argument \code{xNames}
      (in the same order);
      if the elements of \code{coef} have names,
      the element named \code{a_0} is taken as intercept of the \emph{logged}
      equation
      and the elements named \code{a_1}, \ldots, \code{a_n}
      are taken as coefficients of the independent variables
      defined in argument \code{xNames} (numbered in that order).}
   \item{coefCov}{optional covariance matrix of the coefficients
      (the order of the rows and columns must correspond
      to the order of the coefficients in argument \code{coef}).}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   A vector containing the endogenous variable.
   If the inputs are provided as logarithmic values
   (argument \code{dataLogged} is \code{TRUE}),
   the endogenous variable is returned as logarithm;
   non-logarithmic values are returned otherwise.

   If argument \code{coefCov} is specified,
   the returned vector has an attribute \code{"variance"}
   that is a vector containing the variances
   of the calculated (fitted) endogenous variable.
}

\seealso{\code{\link{translogCalc}}, \code{\link{cobbDouglasOpt}}.}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # estimate a Cobb-Douglas production function
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, linear = TRUE )

   # fitted values
   fitted <- cobbDouglasCalc( c( "qLabor", "land", "qVarInput", "time" ), germanFarms,
      coef( estResult )[ 1:5 ] )
   #equal to estResult$fitted

   # fitted values and their variances
   fitted2 <- cobbDouglasCalc( c( "qLabor", "land", "qVarInput", "time" ), germanFarms,
      coef( estResult )[ 1:5 ], coefCov = vcov( estResult )[ 1:5, 1:5 ] )
   # t-values
   c( fitted2 ) / attributes( fitted2 )$variance^0.5
}

\keyword{models}
