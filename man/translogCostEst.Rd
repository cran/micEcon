\name{translogCostEst}
\alias{translogCostEst}
\title{Estimate a translog Cost Function}

\description{
   Estimate a translog cost function.

   NOTE: this function is still under development and incomplete!
}

\usage{
   translogCostEst( cName, yName, pNames, data, fNames = NULL, 
      shifterNames = NULL, dataLogged = FALSE, homPrice = TRUE, \dots )
}

\arguments{
   \item{cName}{a string containing the name of the variable 
      for total cost.}
   \item{yName}{a string containing the name of the variable
      for the total output quantity.}
   \item{pNames}{a vector of strings containing the names
      of the input prices.}
   \item{data}{data frame containing the data
      (possibly a panel data frame created with
      \code{\link[plm]{plm.data}}).}
   \item{fNames}{a vector of strings containing the names
      of fixed inputs.}
   \item{shifterNames}{a vector of strings containing the names
      of the independent variables that should be included
      as shifters only (not in quadratic or interaction terms).}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
   \item{homPrice}{logical. Should homogeneity of degree one in prices be 
      imposed?}
   \item{\dots}{further arguments
      are passed to \code{\link[stats]{lm}}
      or \code{\link[plm]{plm}}.}
}

\value{
   a list of class \code{translogCostEst} containing following objects:
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
   \item{cName}{argument \code{cName}.}
   \item{yName}{argument \code{yName}.}
   \item{pNames}{argument \code{pNames}.}
   \item{fNames}{argument \code{fNames}.}
   \item{shifterNames}{argument \code{shifterNames}.}
   \item{dataLogged}{argument \code{dataLogged}.}
   \item{homPrice}{argument \code{homPrice}.}
}

\seealso{\code{\link{translogEst}}
   and \code{\link{quadFuncEst}}.}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # value of labor input
   germanFarms$vLabor <- germanFarms$pLabor + germanFarms$qLabor
   # total variable cost
   germanFarms$cost <- germanFarms$vLabor + germanFarms$vVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # estimate a translog cost function
   estResult <- translogCostEst( cName = "cost", yName = "qOutput", 
      pNames = c( "pLabor", "pVarInput" ), fNames = "land",
      shifterNames = "time", data = germanFarms, homPrice = FALSE )

   summary( estResult$est )
}

\keyword{models}
