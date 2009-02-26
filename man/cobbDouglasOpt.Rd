\name{cobbDouglasOpt}
\alias{cobbDouglasOpt}
\title{Optimal Values of Independent Variables of a Cobb-Douglas Function}

\description{
   Calculate the optimal values of the variable independent variables
   of a Cobb-Douglas function.
}

\usage{
cobbDouglasOpt( pyName, pxNames, data, coef,
   zNames = NULL, zCoef = NULL, xNames = NULL, dataLogged = FALSE  )
}

\arguments{
   \item{pyName}{character string containing the name of the
      price of the dependent variable.}
   \item{pxNames}{a vector of strings containing the names of the
      prices of the variable independent variables.}
   \item{data}{data frame containing the data.}
   \item{coef}{vector containing the intercept and the coefficients
      of the variable independent variables:
      if the elements of the vector have no names,
      the first element is taken as intercept of the \emph{logged} equation
      and the following elements are taken as coefficients of
      the variable independent variables with corresponding prices
      defined in argument \code{pxNames} (in the same order);
      if the elements of \code{coef} have names,
      the element named \code{a_0} is taken as intercept of the \emph{logged}
      equation
      and the elements named \code{a_1}, \ldots, \code{a_n}
      are taken as coefficients of the variable independent variables
      with corresponding prices
      defined in argument \code{xNames} (numbered in that order).}
   \item{zNames}{optional vector of strings containing the names of
      the fixed independent variables.}
   \item{zCoef}{vector containing the coefficients
      of the fixed independent variables:
      if the elements of the vector have no names,
      they are taken as coefficients of the fixed independent variables
      defined in argument \code{zNames} (in the same order);
      if the elements of \code{coef} have names,
      the elements named \code{d_1}, \ldots, \code{d_m}
      are taken as coefficients of the fixed independent variables
      with corresponding prices
      defined in argument \code{zNames} (numbered in that order).}
   \item{xNames}{optional vector of strings containing the names
      that should be assigned to the returned variable independent variables.}
   \item{dataLogged}{logical. Are the prices and fixed independent variables
      in \code{data}
      with names defined in \code{pyName}, \code{pxNames}, and \code{zNames}
      already logged?}
}

\value{
   A data frame containing the optimal values of the variable
   independent variables.
   If the prices and fixed independent variables are provided
   as logarithmic values
   (argument \code{dataLogged} is \code{TRUE}),
   the optimal values of the variable independent variables are returned
   as logarithms, too;
   non-logarithmic values are returned otherwise.
}

\seealso{\code{\link{cobbDouglasCalc}}.}

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
   estResult <- translogEst( yName = "qOutput",
      xNames = c( "qLabor", "qVarInput", "land", "time" ),
      data = germanFarms, linear = TRUE )

   # calculate optimal quantities of variable inputs
   xCoef <- coef( estResult )[ 1:3 ]
   zCoef <- coef( estResult )[ 4:5 ]
   names( zCoef ) <- c( "d_1", "d_2" )
   optInput <- cobbDouglasOpt( pyName = "pOutput",
      pxNames = c( "pLabor", "pVarInput" ), coef = xCoef,
      data = germanFarms, xNames = c( "qLabor", "qVarInput" ),
      zNames = c( "land", "time" ), zCoef = zCoef )

   # compare observed with optimal input quantities
   plot( germanFarms$qLabor, optInput$qLabor )
   plot( germanFarms$qVarInput, optInput$qVarInput )
}

\keyword{models}
