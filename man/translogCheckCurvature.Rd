\name{translogCheckCurvature}
\alias{translogCheckCurvature}
\alias{print.translogCheckCurvature}
\title{Curvature of a Translog Function}

\description{
   Check curvature of a translog function.
}

\usage{
translogCheckCurvature( xNames, data, coef, convexity = TRUE,
   quasi = FALSE, quadHalf = TRUE, dataLogged = FALSE, ... )

\method{print}{translogCheckCurvature}( x, ... )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients.}
   \item{convexity}{logical. Check whether the function is (quasi)convex
      (default, \code{TRUE}) or (quasi)concave (\code{FALSE}).}
   \item{quasi}{logical. Check whether the function is quasiconvex/quasiconcave
      (\code{TRUE}) or convex/concave (default, \code{FALSE}).}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
   \item{x}{an object returned by \code{translogCheckCurvature}.}
   \item{...}{arguments passed from \code{translogCheckCurvature}
      to \code{\link{semidefiniteness}}
      (if argument \code{quasi} is \code{FALSE}),
      \code{\link{quasiconvexity}}
      (if arguments \code{convexity} and \code{quasi} are both \code{TRUE}),
      or \code{\link{quasiconcavity}}
      (if argument \code{convexity} is \code{FALSE}
      and \code{quasi} is \code{TRUE}).
      Further arguments to \code{print.translogCheckCurvature}
      are currently ignored.}
}

\value{
   \code{translogCheckCurvature} returns a list of class
   \code{translogCheckCurvature} containing following objects:
   \item{obs}{a vector indicating whether the condition
      for the specified curvature is fulfilled at each observation.}
   \item{convexity}{argument \code{convexity}.}
   \item{quasi}{argument \code{quasi}.}
}

\seealso{
   \code{\link{translogEst}} and \code{\link{translogCheckMono}}
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

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

   # check whether the production function is quasiconcave
   translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms, coef( estResult ), convexity = FALSE, quasi = TRUE )
}

\keyword{models}
