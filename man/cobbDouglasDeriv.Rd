\name{cobbDouglasDeriv}
\alias{cobbDouglasDeriv}
\title{Derivatives of a Cobb-Douglas function}

\description{
   Calculate the derivatives of a Cobb-Douglas function.
}

\usage{
cobbDouglasDeriv( xNames, data, coef, coefCov = NULL,
   yName = NULL, dataLogged = FALSE )
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
   \item{yName}{an optional string containing the name of the dependent
      variable.
      If it is \code{NULL}, the dependent variable is calculated
      from the independent variables and the coefficients.}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   a list of class \code{cobbDouglasDeriv} containing following objects:
   \item{deriv}{data frame containing the derivatives.}
   \item{variance}{data frame containing the variances of the derivatives
      (only if argument \code{coefCov} is provided).
      NOTE: if argument \code{yName} is specified,
      the variance of the endogenous variable is currently ignored.}
}

\seealso{\code{\link{cobbDouglasCalc}}, \code{\link{translogDeriv}}.}

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
   estResult <- translogEst( "qOutput", c( "qLabor", "qVarInput", "land", "time" ),
      germanFarms, linear = TRUE )

   # compute the marginal products of the inputs (with "fitted" Output)
   margProducts <- cobbDouglasDeriv( c( "qLabor", "qVarInput", "land", "time" ),
      data = germanFarms, coef = coef( estResult )[1:5],
      coefCov = vcov( estResult )[1:5,1:5] )
   margProducts$deriv
   # t-values
   margProducts$deriv / margProducts$variance^0.5

   # compute the marginal products of the inputs (with observed Output)
   margProductsObs <- cobbDouglasDeriv( c( "qLabor", "qVarInput", "land", "time" ),
      data = germanFarms, coef = coef( estResult )[1:5], yName = "qOutput",
      coefCov = vcov( estResult )[1:5,1:5] )
   margProductsObs$deriv
   # t-values
   margProductsObs$deriv / margProductsObs$variance^0.5
}

\keyword{models}
