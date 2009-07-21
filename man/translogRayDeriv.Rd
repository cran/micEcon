\name{translogRayDeriv}
\alias{translogRayDeriv}
\title{Derivatives of a Translog Ray Function}

\description{
   Calculate the partial derivatives of a translog ray function
   with respect to the independent and the dependent variables.
   Please note that \code{translogRayDeriv} does \emph{not} return
   the partial derivatives of the \sQuote{distance} (\eqn{||y||})
   but of the function \eqn{F(y,x)= Q( \log x, theta ) - \log ||y|| = 0},
   where \eqn{Q(.)} denotes a quadratic function.
}

\usage{
translogRayDeriv( yNames, xNames, data, coef,
   dataLogged = FALSE )
}

\arguments{
   \item{yNames}{a single character string or a vector of character strings
      containing the name(s) of the output quantity.}
   \item{xNames}{a vector of strings containing the names of the
      input quantities.}
   \item{data}{dataframe containing the data.}
   \item{coef}{vector containing all coefficients:
      if there are \code{n} independent variables in \code{xNames},
      the \code{n+2} alpha coefficients must have names
      \code{a_0}, \ldots, \code{a_n}, \code{a_t},
      and the \code{(n+1)*(n+2)/2} beta coefficients must have names
      \code{b_1_1}, \ldots, \code{b_1_n}, \code{b_1_t},
      \code{b_2_2}, \ldots, \code{b_2_n}, \code{b_2_t},
      \ldots, \code{b_n_n}, \code{b_n_t}, \code{b_t_t}
      (only the elements of the upper right triangle of the beta matrix
      are directly obtained from \code{coef};
      the elements of the lower left triangle are obtained by assuming
      symmetry of the beta matrix).}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
}

\value{
   A data frame containing the partial derivatives of the function \eqn{F(y,x)}
   (see above) with respect to \eqn{x} and \eqn{y}.
}

\seealso{
   \code{\link{translogRayEst}}, \code{\link{translogDeriv}}.
}

\author{Arne Henningsen and Geraldine Henningsen}

\examples{
   data( germanFarms )
   # quantity of crop outputs:
   germanFarms$qCrop <- germanFarms$vCrop / germanFarms$pOutput
   # quantity of animal outputs:
   germanFarms$qAnimal <- germanFarms$vAnimal / germanFarms$pOutput
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput

   # estimate a translog ray production function
   estResult <- translogRayEst( yNames = c( "qCrop", "qAnimal" ),
      xNames = c( "qLabor", "land", "qVarInput" ),
      data = germanFarms )

   # compute the partial derivatives
   deriv <- translogRayDeriv( yNames = c( "qCrop", "qAnimal" ),
      xNames = c( "qLabor", "land", "qVarInput" ),
      data = germanFarms, coef = coef( estResult ) )
}

\keyword{models}
