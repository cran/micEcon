\name{translogRayEst}
\alias{translogRayEst}
\alias{print.translogRayEst}
\title{Estimate a Translog Ray Function}

\description{
   Estimate a translog ray function with two endogenous variables,
   e.g.\ a ray production function with two outputs.
}

\usage{
translogRayEst( yNames, xNames, data, shifterNames = NULL, \dots )

\method{print}{translogRayEst}( x, \dots )
}

\arguments{
   \item{yNames}{a vector of two character strings containing
      the names of the two dependent variables.}
   \item{xNames}{a vector of strings containing the names of the
      independent variables.}
   \item{data}{data frame containing the data
      (possibly a panel data frame created with
      \code{\link[plm]{plm.data}}).}
   \item{shifterNames}{a vector of strings containing the names
      of the independent variables that should be included
      as shifters only (not in quadratic or interaction terms).}
   \item{x}{an object of class \code{translogRayEst}.}
   \item{\dots}{further arguments of \code{translogRayEst}
      are passed to \code{\link{translogEst}}
      and possibly further to \code{\link[stats]{lm}}
      or \code{\link[plm]{plm}};
      further arguments of \code{print.translogEst}
      are currently ignored.}
}

\details{
   The actual estimation is done by \code{\link{translogEst}}.
}

\value{
   A list of class \code{translogRayEst} (and \code{translogEst})
   that is very similar to the object returned by \code{\link{translogEst}}.
   It additionally includes following objects:
   \item{distance}{the \dQuote{distance} from the origin (zero)
      to the point of the dependent variables.}
   \item{theta}{the \dQuote{direction} from the origin (zero)
      to the point of the dependent variables.}
}

\seealso{\code{\link{translogEst}}.}

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

   estResult

   summary( estResult )
}

\keyword{models}
