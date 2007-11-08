\name{coef.selection}
\alias{coef.selection}
\alias{coef.summary.selection}
\alias{print.coef.selection}

\title{Extract Coefficients from Selection Models}

\description{
  This function extracts coefficients
  from sample selection models
}

\usage{
\method{coef}{selection}(object, part = "full", ...)
\method{coef}{summary.selection}(object, part = "full", ...)
\method{print}{coef.selection}( x, prefix = TRUE,
      digits = max(3, getOption("digits") - 3), ... )
}

\arguments{
  \item{object}{object of class \code{selection}
    or \code{summary.selection}.}
  \item{part}{character string indication which parts to extract:
   "full" for all estimated parameters (selection estimates,
   outcome estimates, error variance and correlation)
   or "outcome" for the outcome estimates only
   (including the coefficient of the inverse Mill's ratio
   in case of a two-step estimation).}
  \item{x}{object returned by \code{coef.selection}.}
  \item{prefix}{logical. Add a prefix to the names of the coefficients
   that indicates to which equation the coefficient belongs.}
  \item{digits}{numeric, (suggested) number of significant digits.}
  \item{\dots}{currently not used.}
}
\value{
   \code{coef.selection} returns a vector of the estimated coefficients.

   \code{coef.summary.selection} returns a matrix of the estimated coefficients,
      their standard errors, t-values, and p-values.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet (\email{otoomet@ut.ee})}

\seealso{\code{\link{coefficients}}, \code{\link{selection}},
   \code{\link{summary.selection}}, and \code{\link{vcov.selection}}}

\examples{
## Estimate a simple female wage model taking into account the labour
## force participation
   data(Mroz87)
   a <- heckit(lfp ~ huswage + kids5 + mtr + fatheduc + educ + city,
               log(wage) ~ educ + city, data=Mroz87)
## extract all coefficients of the model:
coef( a )

## now extract the coefficients of the outcome model only:
coef( a, part="outcome")

## extract all coefficients, standard errors, t-values
## and p-values of the model:
coef( summary( a ) )

## now extract the coefficients, standard errors, t-values
## and p-values of the outcome model only:
coef( summary( a ), part="outcome")
}

\keyword{methods}
