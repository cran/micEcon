\name{quadFuncEst}
\alias{quadFuncEst}
\title{Estimate a quadratic function}

\description{
   Estimate a quadratic function.
}

\usage{
   quadFuncEst( yName, xNames, data, quadHalf = TRUE, exVarScale = 1 )
}

\arguments{
   \item{yName}{a vector of strings containing the name of the dependent
      variable.}
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
   \item{exVarScale}{a scalar or vector with length equal to nrow( data ).
      All exogenous variables (including the quadratic terms) are
      divided by this value/these values.}
}

\value{
   a list of class \code{quadFuncEst} containing following objects:
   \item{lm}{a list returned by \code{lm}.}
   \item{residuals}{residuals.}
   \item{fitted}{fitted values.}
   \item{coef}{a list containing the the estimated coefficients:
      alpha0 = \eqn{\alpha_0} (scalar),
      alpha = \eqn{\alpha_i} (vector),
      beta = \eqn{\beta_ij} (matrix).}
   \item{allCoef}{vector of all coefficients.}
   \item{allCoefCov}{Covariance matrix of all coefficients.}
   \item{r2}{\eqn{R^2} value.}
   \item{r2bar}{adjusted \eqn{R^2} value.}
   \item{nObs}{number of observations.}
   \item{model.matrix}{the model matrix.}
}

\seealso{\code{\link{quadFuncCalc}}, \code{\link{quadFuncDeriv}},
   \code{\link{translogEst}} and \code{\link{snqProfitEst}}.}

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

   estResult$coef
   estResult$r2
}

\keyword{models}
