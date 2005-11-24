\name{aidsTestConsist}
\alias{aidsTestConsist}
\title{Consistency Test of the AIDS}
\description{
   Test whether the specification of the AIDS is consistent
   with microeconomic demand theory (i.e. utility maximisation).
}

\usage{
   aidsTestConsist( pNames, wNames, xtName, data, coef = NULL,
      alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) )
}

\arguments{
   \item{pNames}{a vector of strings containing the names of the prices.}
   \item{wNames}{a vector of strings containing the names of the expenditure
      shares.}
   \item{xtName}{a string containing the variable name of total expenditure.}
   \item{data}{a data frame containing the data.}
   \item{coef}{a list containing the coefficients alpha, beta and gamma.}
   \item{alpha0}{coefficient \eqn{\alpha_0} of the translog price index.}
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

\seealso{\code{\link{aidsEst}}, \code{\link{aidsEla}}}

\examples{
   data( Blanciforti86 )
   # Data on food consumption are available only for the first 32 years
   Blanciforti86 <- Blanciforti86[ 1:32, ]

   pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   estResult <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86, method = "IL:L" )
   tc <- aidsTestConsist( pNames, wNames, "xFood", Blanciforti86,
      coef = estResult$coef )
   tc$mPercent   # great!
   tc$cPercent   # Oh, that is bad!
}

\keyword{models}
