\name{aidsBestA0}
\alias{aidsBestA0}
\title{Find 'best' Value for Coefficient alpha 0}

\description{
   Search for the intercept of the translog price index
   (\eqn{\alpha_0}) that gives the best fit to the model,
   i.e. the value that leads to the smallest determinant
   of the residual covariance matrix (see Michalek
   and Keyzer, 1992).
}

\usage{
   aidsBestA0( pNames, wNames, xtName,
      data = NULL, ivNames = NULL, method = "MK:L",
      a0min = -50, a0max = 50, stoprange = 3, stopiter = 10,
      verbose = FALSE, ... )
}

\arguments{
   \item{pNames}{a vector of strings containing the names of prices.}
   \item{wNames}{a vector of strings containing the names of expenditure
      shares.}
   \item{xtName}{a name of the total expenditure variable.}
   \item{data}{a data frame containing the data.}
   \item{ivNames}{a vector of strings containing the names of instrumental
      variables.}
   \item{method}{the method to estimate the aids (see \code{\link{aidsEst}}).}
   \item{a0min}{lower limit of the range for searching for \eqn{\alpha_0}.}
   \item{a0max}{upper limit of the range for searching for \eqn{\alpha_0}.}
   \item{stoprange}{stop searching when the search interval is smaller than
      or equal to \code{stoprange}.}
   \item{stopiter}{maximal number of iterations.}
   \item{verbose}{print each determinant of the residual covariance matrix
      immediately after its calculation.}
   \item{...}{argunents passed to \code{\link{aidsEst}}.}
}

\details{
   Since the call of \code{\link{aidsEst}} takes some time,
   the search algorithm is constructed to minimize the
   calls of the function \code{\link{aidsEst}}.
}

\value{
   a list containing following objects:
  \item{alpha0}{\eqn{\alpha_0} that gives the best fit.}
  \item{allValues}{all \eqn{\alpha_0} tested and the corresponding
      determinants of the residual covariance matrix.}
  \item{iter}{number of iterations.}
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Michalek, J. and M. A. Keyzer (1992)
   Estimation of a two-stage LES-AIDS consumer demand system
   for eight EC countries.
   \emph{European Review of Agricultural Economics}, 19 (2), p. 137-163.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{ \code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   bestA0 <- aidsBestA0( c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), "xFood",
      data = Blanciforti86, method = "MK:L" ) # may take some time
   bestA0$alpha0
   plot( bestA0$allValues ) # this should be convex
}

\keyword{ models }
