\name{translogEst}
\alias{translogEst}
\title{Estimate a translog function}

\description{
   Estimate a translog function.
}

\usage{
   translogEst( yName, xNames, data, quadHalf = TRUE, logValues = FALSE )
}

\arguments{
   \item{yName}{a string containing the name of the dependent
      variable.}
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
   \item{logValues}{logical. Are the values in \code{data} already logged?}
}

\value{
   a list of class \code{translogEst} containing following objects:
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

\seealso{\code{\link{translogCalc}}, \code{\link{translogDeriv}}
   and \code{\link{quadFuncEst}}.}

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

   estResult$coef
   estResult$r2
}

\keyword{models}
