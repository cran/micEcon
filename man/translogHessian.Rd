\name{translogHessian}
\alias{translogHessian}
\title{Hessian matrix of a translog function}

\description{
   Calculate the Hessian matrices of a translog function.
}

\usage{
   translogHessian( xNames, data, coef, yName = NULL,
   dataLogged = FALSE, bordered = FALSE )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients.}
   \item{yName}{an optional string containing the name of the dependent
      variable.
      If it is \code{NULL}, the dependent variable is calculated
      from the independent variables and the coefficients.}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
   \item{bordered}{logical. Should the \emph{bordered} Hessians be returned?}
}

\value{
   a list containing following the (bordered) Hessian matrices at
   each data point.
}

\seealso{
   \code{\link{translogEst}}, \code{\link{translogDeriv}} and
   \code{\link{translogCalc}}
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

   # estimate a quadratic production function
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms )

   # compute the Hessian matrices (with "fitted" output)
   hessians <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, coef( estResult ) )
   hessians[[ 1 ]]

   # compute the Hessian matrices (with observed output)
   hessiansObs <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, coef( estResult ), "qOutput" )
   hessiansObs[[ 1 ]]

   # compute the bordered Hessian matrices
   borderedHessians <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, coef( estResult ), bordered = TRUE )
   borderedHessians[[ 1 ]]
}

\keyword{models}
