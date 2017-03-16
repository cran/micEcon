\name{quadFuncEst}
\alias{quadFuncEst}
\title{Estimate a quadratic function}

\description{
   Estimate a quadratic function.
}

\usage{
quadFuncEst( yName, xNames, data, shifterNames = NULL,
   linear = FALSE, homWeights = NULL, 
   regScale = 1, \dots )
}

\arguments{
   \item{yName}{a character string containing the name of the dependent
      variable.}
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{data frame containing the data
      (possibly a panel data frame created with
      \code{\link[plm]{pdata.frame}}).}
   \item{shifterNames}{a vector of strings containing the names
      of the independent variables that should be included
      as shifters only (not in quadratic or interaction terms).}
   \item{linear}{logical. Restrict the coefficients of all quadratic
      and interaction terms to be zero so that the estimated function
      is linear in the exogenous variables?}
   \item{homWeights}{numeric vector with named elements that are weighting factors
      for calculating an index that is used to normalize the variables
      for imposing homogeneity of degree zero in these variables
      (see details).}
   \item{regScale}{a scalar or vector with length equal to \code{nrow( data )}.
      All regressors except for shifter variables that are logical 
      or factors are divided by \code{regScale}
      (NOTE: quadratic and interaction terms are also divided by 
      \code{regScale} and NOT divided by the square of \code{regScale}).}
   \item{\dots}{further arguments are passed to \code{\link[stats]{lm}}
      or \code{\link[plm]{plm}}.}
}

\details{
   If argument \code{homWeights} is used to impose homogeneity of degree
   zero in some variables, 
   the weighting factors in this vector must have names 
   that are equal to the variable names in argument \code{xNames}.
   The order of the elements in \code{homWeights}
   is arbitrary and may or may not be equal 
   to the order of the elements in \code{xNames}.
   Argument \code{homWeights} may contain less elements 
   than \code{xNames};
   in this case, homogeneity of degree zero is imposed only
   on variables with names in \code{homWeights}.
   Please note that the weighting factor of a variable
   (\eqn{P_i}) in \code{homWeights}
   (\eqn{w_i = \partial P / \partial P_i}{w_i = d P / d P_i})
   is not really its weight
   (\eqn{( \partial P / \partial P_i ) ( P_i / P )}{%
      ( d P / d P_i ) ( P_i / P )}),
   in particular,
   if the numerical values of the variables (\eqn{P_1, \ldots, P_n})
   are rather different.
}

\value{
   a list of class \code{quadFuncEst} containing following objects:
   \item{est}{the object returned by \code{\link[stats]{lm}}
      or \code{\link[plm]{plm}}.}
   \item{nExog}{length of argument \code{xNames}.}
   \item{nShifter}{length of argument \code{shifterNames}.}
   \item{residuals}{residuals.}
   \item{fitted}{fitted values.}
   \item{coef}{vector of all coefficients.}
   \item{coefCov}{covariance matrix of all coefficients.}
   \item{r2}{\eqn{R^2} value.}
   \item{r2bar}{adjusted \eqn{R^2} value.}
   \item{nObs}{number of observations.}
   \item{model.matrix}{the model matrix.}
   \item{call}{the matched call.}
   \item{yName}{argument \code{yName}.}
   \item{xNames}{argument \code{xNames}.}
   \item{shifterNames}{argument \code{shifterNames}.}
   \item{homWeights}{argument \code{homWeights}.}
   \item{regScale}{argument \code{regScale}.}
}

\seealso{\code{\link{quadFuncCalc}}, \code{\link{quadFuncDeriv}},
   \code{\link{translogEst}} and \code{\link[micEconSNQP]{snqProfitEst}}.}

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

   coef( estResult )
   estResult$r2
}

\keyword{models}
