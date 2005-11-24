\name{aidsPx}
\alias{aidsPx}
\title{ Price Index for the AIDS }

\description{
   Calculate log price indices for the AIDS and LA-AIDS.
}

\usage{
   aidsPx( px, pNames, wNames = NULL, data = NULL, base = 1, coef = NULL,
      alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) )
}

\arguments{
   \item{px}{the prices index to be used (see details).}
   \item{pNames}{a vector of strings containing the names of the prices.}
   \item{wNames}{a vector of strings containing the names of the expenditure
      shares.}
   \item{data}{a data frame containing the data.}
   \item{base}{the base to calculate the indices (see details).}
   \item{coef}{a list containing the coefficients alpha, beta and gamma
      (only needed for the translog price index).}
   \item{alpha0}{coefficient \eqn{\alpha_0} of the translog price index.}
}
\details{
   Currently 6 different price indices are implemented:
   \begin{itemize}
   \item translog price index ('TL'),
   \item Stone price index ('S'),
   \item Stone price index with lagged shares ('SL'),
   \item loglinear analogue to the Paasche price index ('P'),
   \item loglinear analogue of the Laspeyres price index ('L'), and
   \item Tornqvist price index ('T').
   \end{itemize}

   The argument \code{base} can be either
   \begin{itemize}
   \item a single number: the row number of the base prices and quantities,
   \item a vector indicating several observations: The means of these
      observations are used as base prices and quantities, or
   \item a logical vector with the same length as the \code{data}: The
      means of the observations indicated as 'TRUE' are used as base prices
      and quantities.
   \end{itemize}
}

\value{
   a vector containing the log price index.
}

\references{
   Deaton, A.S. and J. Muellbauer (1980)
   An Almost Ideal Demand System.
   \emph{American Economic Review}, 70, p. 312-326.

   Moschini, G. (1995)
   Units of Measurement and the Stone Index in Demand System Estimation.
   \emph{American Journal of Agricultural Economics}, 77, p. 63-68.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{ \code{\link{aidsEst}}}

\examples{
   data( Blanciforti86 )
   aidsPx( "S", c( "pFood1", "pFood2", "pFood3", "pFood4" ),
      c( "wFood1", "wFood2", "wFood3", "wFood4" ), Blanciforti86 )
}

\keyword{ models }
