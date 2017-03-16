\name{translogProdFuncMargCost}
\alias{translogProdFuncMargCost}
\title{Marginal Costs of Translog Production Function}

\description{
   Calculate the marginal costs of the output
   from a translog production function.
}

\usage{
translogProdFuncMargCost( yName, xNames, wNames, data, coef,
   dataLogged = FALSE )
}

\arguments{
   \item{yName}{a single character string
      containing the name of the output quantity.}
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
   A vector containing the marginal costs of producing the output.
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
   margCost <- translogProdFuncMargCost( yName = "qOutput",
      xNames = c( "qLabor", "land", "qVarInput", "time" ),
      wNames = c( "pLabor", "pLand", "pVarInput", "pTime" ),
      data = germanFarms, coef = coef( estResult ) )
}

\keyword{models}
