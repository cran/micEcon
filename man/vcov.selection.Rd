\name{vcov.selection}
\alias{vcov.selection}

\title{Extract Variance Covariance Matrix}
\description{
   This function extracts the coefficient variance-covariance matrix
   from sample selection models.
}
\usage{
\method{vcov}{selection}(object, part = "full", ...)
}
\arguments{
   \item{object}{object of class "selection".}
   \item{part}{character string indication which parts
      of the variance-covariance matrix to extract:
      "full" for all estimated parameters (selection estimates,
      outcome estimates, error variance and correlation)
      or "outcome" for the outcome estimates only
      (including the coefficient of the inverse Mill's ratio
      in case of a two-step estimation).}
   \item{\dots}{currently not used.}
}
\details{
   The variance-covariance matrix of a two-step estimate
   is currently only partly implemented.
   The unimplemented part of the matrix is filled with NAs.
}
\value{
   the estimated variance covariance matrix of the coefficients.
}
\author{
   Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet \email{otoomet@ut.ee}
}

\seealso{\code{\link[stat]{vcov}}, \code{\link{selection}},
   \code{\link{coef.selection}}.}
\examples{
## Estimate a simple female wage model taking into account the labour
## force participation
   data(Mroz87)
   a <- heckit(lfp ~ huswage + kids5 + mtr + fatheduc + educ + city,
               log(wage) ~ educ + city, data=Mroz87)
## extract the full variance-covariance matrix:
vcov( a )
## now extract the variance-covariance matrix of the outcome model only:
vcov( a, part = "outcome" )
}
\keyword{methods}
