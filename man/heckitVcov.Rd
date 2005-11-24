\name{heckitVcov}
\alias{heckitVcov}
\title{Heckit Variance Covariance Matrix}

\description{
   Calculate Variance Covariance Matrix for Heckit Estimations.
}

\usage{
   heckitVcov( xMat, wMat, vcovProbit, rho, delta, sigma,
   saveMemory = TRUE )
}

\arguments{
   \item{xMat}{model matrix of the 2nd step estimation.}
   \item{wMat}{model matrix of the 1st step probit estimation.}
   \item{vcovProbit}{variance covariance matrix of the 1st step
      probit estimation.}
   \item{rho}{the estimated \eqn{\rho}, see Greene (2003, p. 784).}
   \item{delta}{the estimated \eqn{\delta}s, see Greene (2003, p. 784).}
   \item{sigma}{the estimated \eqn{\sigma}, see Greene (2003, p. 784).}
   \item{saveMemory}{logical. Save memory by using a different implementation
      of the formula? (this should not influence the results).}
}

\value{
   the variance covariance matrix of the coefficients.
}

\references{
  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Lee, L., G. Maddala and R. Trost (1980)
  Asymetric covariance matrices of two-stage probit and two-stage tobit
  methods for simultaneous equations models with selectivity.
  \emph{Econometrica}, 48, p. 491-503.
}

\seealso{\code{\link{heckit}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\keyword{models}
