\name{aidsTestConsist}
\alias{aidsTestConsist}
\alias{testConsist.aidsEst}
\title{Consistency Test of the AIDS}
\description{
   Test whether the specification of the AIDS is consistent
   with microeconomic demand theory (i.e. utility maximisation).
}

\usage{
   aidsTestConsist( priceNames, shareNames, totExpName, data, coef = NULL,
      alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) )

   \method{testConsist}{aidsEst}( object, ... )
}

\arguments{
   \item{priceNames}{a vector of strings containing the names of the prices.}
   \item{shareNames}{a vector of strings containing the names of the expenditure
      shares.}
   \item{totExpName}{a string containing the variable name of total expenditure.}
   \item{data}{a data frame containing the data.}
   \item{coef}{a list containing the coefficients alpha, beta and gamma.}
   \item{alpha0}{coefficient \eqn{\alpha_0} of the translog price index.}
   \item{object}{an object of class \code{aidsEst}.}
   \item{...}{currently not used.}
}

\details{
   \code{testConsist.aidsEst} is a wrapper function to \code{aidsTestConsist}
   that extracts all relevant arguments for \code{aidsTestConsist}
   from \code{object}.
}

\value{
   a list containing following elements:
   \item{monotony}{a logical vector indicating whether the monotony
      condition is fulfilled at each observation.}
   \item{mPercent}{percent of observations where the monotony condition is
      fulfilled.}
   \item{concavity}{a logical vector indicating whether the concavity
      condition is fulfilled at each observation.}
   \item{cPercent}{percent of observations where the concavity condition is
      fulfilled.}
   \item{cMatrices}{a list of the 'C' matrices for each observation to check
      for concavity (see Deaton and Muellbauer, 1980b, p.76 ).}
}
\references{
   Deaton, A.S. and J. Muellbauer (1980a)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Deaton, A.S. and J. Muellbauer (1980b)
   \emph{Economics and Consumer Behavior},
   Cambridge University Press, Cambridge.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{aidsEst}}, \code{\link{aidsElas}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   priceNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   shareNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   estResult <- aidsEst( priceNames, shareNames, "xFood",
      data = Blanciforti86, method = "IL:L" )
   tc <- aidsTestConsist( priceNames, shareNames, "xFood", Blanciforti86,
      coef = estResult$coef )
   tc$mPercent   # great!
   tc$cPercent   # Oh, that is bad!

   # the same can be obtained in an easier way
   testConsist.aidsEst( estResult )
}

\keyword{models}
