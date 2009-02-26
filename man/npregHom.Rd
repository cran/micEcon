\name{npregHom}
\alias{npregHom}
\title{Nonparametric Regression with Homogeneity Imposed}

\description{
   Nonparametric regression with homogeneity of degree zero
   in some regressors imposed.
}

\usage{
npregHom( yName, xNames, homWeights, data, restrictGrad = TRUE,
   bws = NULL, \dots )
}

\arguments{
   \item{yName}{a character string containing the name of the dependent
      variable.}
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{homWeights}{numeric vector with named elements that are weighting factors
      for calculating an index that is used to normalize the variables
      for imposing homogeneity of degree zero in these variables
      (see details).}
   \item{data}{data frame containing the data.}
   \item{restrictGrad}{logical value indicating
      whether the sum of the gradients of all normalized variables
      should be restricted to be zero? (see details).}
   \item{bws}{bandwidths (see \code{\link[np]{npreg}}).}
   \item{\dots}{further arguments are passed to \code{\link[np]{npreg}}.}
}

\details{
   Argument \code{homWeights} is used to impose homogeneity of degree
   zero in some (continuous) independent variables.
   The weighting factors in this vector must have names
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

   The variables that are normalized with the weighted index of theses variables
   are linearly dependent.
   Hence, a model that includes these variables cannot be estimated
   by standard econometric methods such as OLS.
   To allow the estimation of this model by standard econometric methods,
   the sum of the gradients (=coefficients) of the normalized variables
   is generally restricted to zero.
   If argument \code{restrictGrad} is \code{TRUE},
   this is done also by \code{npregHom}.
   In contrast to OLS results that do not depend on
   which variable is eliminated by the restriction,
   the results of \code{npregHom} depend on
   which variable is eliminated by the restriction.
   The variable that corresponds to the first weight
   in argument \code{homWeights} is eliminated in \code{npregHom}.
}

\value{
   a list of class \code{npregHom} containing following objects:
   \item{est}{the object returned by \code{\link[np]{npreg}}.}
   \item{grad}{matrix containing the gradients of all regressors.}
   \item{call}{the matched call.}
   \item{yName}{argument \code{yName}.}
   \item{xNames}{argument \code{xcNames}.}
   \item{homWeights}{argument \code{homWeights}.}
   \item{restrictGrad}{argument \code{restrictGrad}.}
}

\seealso{\code{\link{elas.npregHom}}, \code{\link{quadFuncEst}},
   \code{\link[np]{npreg}}, and \code{\link[np]{npregbw}}.}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # weights to impose
   weights <- c(
      pOutput = mean( germanFarms$qOutput ),
      pVarInput = mean( germanFarms$qVarInput ),
      pLabor = mean( germanFarms$qLabor ) )
   weights <- weights / sum( weights )

   # estimate an input demand function
   estResult <- npregHom( "qVarInput",
      xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
      data = germanFarms, homWeights = weights )
   estResult$grad

   # estimate an input demand function using the Epanechnikov kernel
   estResultEpa <- npregHom( "qVarInput",
      xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
      data = germanFarms, homWeights = weights, ckertype="epanechnikov" )
   estResultEpa$grad

   # estimate an input demand function with manual bandwidths selection
   estResultMan <- npregHom( "qVarInput",
      xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
      data = germanFarms, homWeights = weights, bws = rep( 1, 3 ),
      bwscaling = TRUE )
   estResultMan$grad
}

\keyword{models}
\keyword{nonparametric}
\keyword{regression}
