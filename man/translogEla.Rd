\name{translogEla}
\alias{translogEla}
\alias{elas.translogEst}
\title{Elasticities of a translog Function}

\description{
   Calculate the elasticities of a translog function.
}

\usage{
translogEla( xNames, data, coef, coefCov = NULL,
   dataLogged = FALSE  )

\method{elas}{translogEst}( object, data = NULL, dataLogged = NULL,
   \dots )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data;
      if argument \code{data} of \code{elas.translogEst} is not specified,
      the data frame that was used for the estimation 
      is used for calculating elasticities.}
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
   \item{dataLogged}{logical. Are the values in \code{data} already logged?
      If argument \code{dataLogged} of \code{elas.translogEst} is not specified,
      the same value as used in \code{translogEst} for creating \code{object}
      is used.}
   \item{object}{object of class \code{translogEst} 
      (returned by \code{\link{translogEst}}).}
   \item{\ldots}{currently ignored.}
}

\details{
   Shifter variables do not need to be specified, 
   because they have no effect on the elasticities.
   Hence, you can use this function to calculate elasticities
   even for translog functions that have been estimated
   with shifter variables.
}

\value{
   A data frame containing the elasticities,
   where each column corresponds to one of the independent variables.
   If argument \code{coefCov} is provided, it has the attributes
   \code{variance} and \code{stdDev},
   which are two data frames containing the variances
   and the standard deviations, respectively, of the elasticities.
}

\seealso{\code{\link{translogEst}} and \code{\link{translogCalc}}}

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

   # calculate production elasticities of all inputs
   estEla <- translogEla( c( "qLabor", "land", "qVarInput", "time" ), 
      data = germanFarms, coef = coef( estResult ), 
      coefCov = vcov( estResult ) )
   # all elasticities
   estEla
   # t-values of all elasticities
   estEla / attributes( estEla )$stdDev
}

\keyword{models}
