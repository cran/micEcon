\name{aidsCalc}
\alias{aidsCalc}
\title{ Shares and Quantities of the Almost Ideal Demand System }

\description{
  Given prices, total expenditure and coefficients this function
  calculates the demanded expenditure shares and quantities based
  on the Almost Ideal Demand System.
}

\usage{
   aidsCalc( pNames, xtName, data = NULL, px = "TL", lnp = NULL,
      coef = NULL, alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) )
}

\arguments{
   \item{pNames}{vector containing the names of the prices. }
   \item{xtName}{name of the total expenditure variable. }
   \item{data}{data frame containing the variables. }
   \item{px}{prices index (see \code{\link{aidsPx}}) }
   \item{lnp}{log of the price index.}
   \item{coef}{a list containing the coefficients alpha, beta and gamma.}
   \item{alpha0}{coefficient \eqn{\alpha_0} in the translog price index.}
}

\value{
   \code{aidsCalc} returns a list with following elements:
   \item{shares}{a data frame containing the calculated expenditure shares}
   \item{quantities}{a data frame containing the calculated quantites}
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{ \code{\link{aidsEst}}, \code{\link{aidsPx}} }

\examples{
   data( Blanciforti86 )
   pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
   wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )

   ## LA-AIDS
   estResult <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86, method = "LA:L" )

   lnp <- aidsPx( "L", pNames, wNames, Blanciforti86 )

   fitted <- aidsCalc( pNames, "xFood", Blanciforti86,
      coef = estResult$coef, lnp = lnp )

   fitted$shares  # equal to estResult$wFitted
   fitted$quant   # equal to estResult$qFitted

   ## AIDS
   estResult <- aidsEst( pNames, wNames, "xFood",
      data = Blanciforti86, method = "MK:L" )

   fitted <- aidsCalc( pNames, "xFood", Blanciforti86,
      coef = estResult$coef )

   fitted$shares  # equal to estResult$wFitted
   fitted$quant   # equal to estResult$qFitted

}

\keyword{ models }
