\name{translogDeriv}
\alias{translogDeriv}
\title{Derivatives of a translog function}

\description{
   Calculate the derivatives of a translog function.
}

\usage{
translogDeriv( xNames, data, coef, coefCov = NULL,
   yName = NULL, dataLogged = FALSE )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients.}
   \item{coefCov}{optional covariance matrix of the coefficients.}
   \item{yName}{an optional string containing the name of the dependent
      variable.
      If it is \code{NULL}, the dependent variable is calculated
      from the independent variables and the coefficients.}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   a list of class \code{translogDeriv} containing following objects:
   \item{deriv}{data frame containing the derivatives.}
   \item{variance}{data frame containing the variances of the derivatives
      (not implemented yet).}
   \item{stdDev}{data frame containing the standard deviations of the derivatives
      (not implemented yet).}
}

\seealso{
   \code{\link{translogEst}}, \code{\link{translogCalc}} and
   \code{\link{translogHessian}}
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

   # estimate a translog production function
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms )

   # compute the marginal products of the inputs (with "fitted" Output)
   margProducts <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, coef( estResult ), vcov( estResult ) )
   margProducts$deriv

   # compute the marginal products of the inputs (with observed Output)
   margProductsObs <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, coef( estResult ), vcov( estResult ), "qOutput" )
   margProductsObs$deriv
}

\keyword{models}
