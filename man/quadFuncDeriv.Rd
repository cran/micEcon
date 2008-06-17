\name{quadFuncDeriv}
\alias{quadFuncDeriv}
\title{Derivatives of a quadratic function}

\description{
   Calculate the derivatives of a quadratic function.
}

\usage{
   quadFuncDeriv( xNames, data, coef, coefCov = NULL, quadHalf = TRUE  )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients.}
   \item{coefCov}{optional covariance matrix of the coefficients.}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
}

\value{
   a list of class \code{quadFuncDeriv} containing following objects:
   \item{deriv}{data frame containing the derivatives.}
   \item{variance}{data frame containing the variances of the derivatives
      (only if argument \code{coefCov} is provided).}
   \item{stdDev}{data frame containing the standard deviations of the derivatives
      (only if argument \code{coefCov} is provided).}
}

\seealso{\code{\link{quadFuncEst}} and \code{\link{quadFuncCalc}}}

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

   # compute the marginal products of the inputs
   margProducts <- quadFuncDeriv( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, coef( estResult ), vcov( estResult ) )
   margProducts$deriv
}

\keyword{models}
