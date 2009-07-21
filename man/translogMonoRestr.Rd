\name{translogMonoRestr}
\alias{translogMonoRestr}

\title{Monotonicity Restrictions of a Translog Function}

\description{
   Create matrix to check or impose the monotonicity restrictions 
   of a translog function.
}

\usage{
translogMonoRestr( xNames, data, quadHalf = TRUE,
   dataLogged = FALSE, box = FALSE )
}

\arguments{
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{dataframe containing the data.}
   \item{quadHalf}{logical. Multiply the quadratic terms by one half?}
   \item{dataLogged}{logical. Are the values in \code{data} already logged?}
   \item{box}{logical. Should monotonicity be imposed within 
      an \eqn{n}-dimensional box that includes all points in \code{data}?
      If \code{FALSE}, monotonicity is imposed (only) within 
      an \eqn{n}-dimensional polygon that includes all points in \code{data}.
      (\eqn{n} is the number of independent varables.)}
}

\value{
   \code{translogMonoRestr} returns a matrix of dimension
   \eqn{( n \cdot N ) \times c},
   where \eqn{n} is the number of independent varables,
   \eqn{N} is the number of data points at which monotonicity
   should be imposed (if argument \code{box} is \code{FALSE},
   \eqn{N} is the number of rows in \code{data};
   if argument \code{box} is \code{TRUE}, \eqn{N = 2^n}),
   and \eqn{c = 1 + n ( n + 3 ) / 2} is the number of (linearly independent) 
   coefficients.
   Multiplying a row of this matrix (e.g.\ the \eqn{k}th row of \eqn{M})
   by the vector of coefficients (\eqn{\beta})
   results in the derivative of the dependent variable (\eqn{y}) 
   with respect to one independent variable (e.g.\ \eqn{x_i})
   at one data point (e.g.\ \eqn{j}): 
   \deqn{M[k,] \cdot \beta = \frac{ \partial \ln y }{ \partial \ln x_i}}{
   M[k,] * \beta = ( d log y )/( d log x_i )},
   evaluated at \eqn{x_{1j}}, \ldots, \eqn{x_{nj}},
   where \eqn{k = ( i - 1 ) N + j}.
   Hence, the observations run faster than the independent variables.
}

\seealso{
   \code{\link{translogEst}}, \code{\link{translogDeriv}},
   and \code{\link{translogCheckMono}}
}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   # quantity of variable inputs
   germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput

   # matrix to check or impose monotonicity at all observations
   monoRestrObs <- translogMonoRestr( c( "qLabor", "land", "qVarInput" ),
      germanFarms )

   # matrix to check or impose monotonicity within a box that includes all 
   # observations
   monoRestrBox <- translogMonoRestr( c( "qLabor", "land", "qVarInput" ),
      germanFarms, box = TRUE )
}

\keyword{models}
