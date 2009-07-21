\name{translogProdFuncMargCost}
\alias{translogProdFuncMargCost}
\title{Marginal Costs of Translog Production Function}

\description{
   Calculate the marginal costs of the output(s)
   from a single-output translog production function or
   a multiple-output translog ray production function.
}

\usage{
translogProdFuncMargCost( yNames, xNames, wNames, data, coef,
   dataLogged = FALSE )
}

\arguments{
   \item{yNames}{a single character string or a vector of character strings
      containing the name(s) of the output quantity.}
   \item{xNames}{a vector of strings containing the names of the
      input quantities.}
   \item{wNames}{a vector of strings containing the names of the
      input prices.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients:
      if there are \code{n} inputs in \code{xNames},
      the \code{n+1} alpha coefficients must have names
      \code{a_0}, \ldots, \code{a_n}
      and the \code{n*(n+1)/2} beta coefficients must have names
      \code{b_1_1}, \ldots, \code{b_1_n}, \ldots, \code{b_n_n}
      (only the elements of the upper right triangle of the beta matrix
      are directly obtained from \code{coef};
      the elements of the lower left triangle are obtained by assuming
      symmetry of the beta matrix).}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   A data frame containing the marginal costs of producing the output.
}

\seealso{
   \code{\link{translogEst}}, \code{\link{translogCalc}},
   \code{\link{translogDeriv}}, \code{\link{translogEla}} and
   \code{\link{translogCostEst}}.
}

\author{Arne Henningsen and Geraldine Henningsen}

\examples{
   data( germanFarms )
   # output quantity:
   germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
   # quantity of crop outputs:
   germanFarms$qCrop <- germanFarms$vCrop / germanFarms$pOutput
   # quantity of animal outputs:
   germanFarms$qAnimal <- germanFarms$vAnimal / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
   # a time trend to account for technical progress:
   germanFarms$time <- c(1:20)

   # generate (artificial) prices
   germanFarms$pLand <- 200 + 15 * germanFarms$time
   germanFarms$pTime <- 1

   # estimate a single-output translog production function
   estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
      germanFarms )

   # compute the marginal costs of producing the output
   margCost <- translogProdFuncMargCost( yNames = "qOutput",
      xNames = c( "qLabor", "land", "qVarInput", "time" ),
      wNames = c( "pLabor", "pLand", "pVarInput", "pTime" ),
      data = germanFarms, coef = coef( estResult ) )

   # estimate a two-output translog ray production function
   estResultRay <- translogRayEst( yNames = c( "qCrop", "qAnimal" ),
      xNames = c( "qLabor", "land", "qVarInput" ),
      data = germanFarms )

   # compute the marginal costs of producing the two outputs
   margCostRay <- translogProdFuncMargCost( yNames = c( "qCrop", "qAnimal" ),
      xNames = c( "qLabor", "land", "qVarInput" ),
      wNames = c( "pLabor", "pLand", "pVarInput" ),
      data = germanFarms, coef = coef( estResultRay ) )
}

\keyword{models}
