\name{quadFuncDeriv}
\alias{quadFuncDeriv}
\title{Derivatives of a quadratic function}

\description{
   Calculate the derivatives of a quadratic function.
}

\usage{
   quadFuncDeriv( xNames, data, coef, coefCov = NULL, 
      homWeights = NULL, quadHalf = TRUE  )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe or a vector with named elements containing the data.}
   \item{coef}{vector containing all coefficients:
      if there are \code{n} exogenous variables in \code{xNames},
      the \code{n+1} alpha coefficients must have names 
      \code{a_0}, \ldots, \code{a_n}
      and the \code{n*(n+1)/2} beta coefficients must have names
      \code{b_1_1}, \ldots, \code{b_1_n}, \ldots, \code{b_n_n}
      (only the elements of the upper right triangle of the beta matrix
      are directly obtained from \code{coef};
      the elements of the lower left triangle are obtained by assuming
      symmetry of the beta matrix).}
   \item{coefCov}{optional covariance matrix of the coefficients:
      the row names and column names must be the same as the names
      of \code{coef}.}
   \item{homWeights}{numeric vector with named elements that are weighting factors
      for calculating an index that is used to normalize the variables
      for imposing homogeneity of degree zero in these variables
      (see documentation of \code{\link{quadFuncEst}}).}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
}

\details{
   Shifter variables do not need to be specified, 
   because they have no effect on the partial derivatives.
   Hence, you can use this function to calculate partial derivatives
   even for quadratic functions that have been estimated
   with shifter variables.
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

\author{Arne Henningsen}

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
