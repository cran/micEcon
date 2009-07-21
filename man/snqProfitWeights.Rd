\name{snqProfitWeights}
\alias{snqProfitWeights}
\title{SNQ Profit function: Weights of prices for normalization}

\description{
   Returns a vector of weights to normalize prices on a
   Symmetric Normalized  Quadratic (SNQ) Profit function.
}

\usage{snqProfitWeights( priceNames, quantNames, data, method = "DW92", base = 1 )}

\arguments{
   \item{priceNames}{a vector of strings containing the names of netput prices.}
   \item{quantNames}{a vector of strings containing the names of netput quantities.}
   \item{data}{a data frame containing the data.}
   \item{method}{the method to determine the weights (see details).}
   \item{base}{the base period(s) for scaling prices (see details).}
}

\details{
   If argument \code{method} is 'DW92' the method of Diewert and
   Wales (1992) is applied. They predetermine the weights by
   \deqn{ \theta_{i} = \frac{
   \displaystyle \left | \overline{x}_{i} \right| p_{i}^{0} }{
   \displaystyle \sum_{i=1}^{n} \left| \overline{x}_{i} \right| p_{i}^{0}}}

   Defining the scaled netput quantities as
   \eqn{\widetilde{x}_{i}^{t} = x_{i}^{t}\cdot p_{i}^{0}}
   we get following formula:

   \deqn{ \theta_{i} = \frac{
   \displaystyle  \left| \overline{ \widetilde{ x } }_{i} \right|}{
   \displaystyle \sum_{i=1}^{n} \left| \overline{ \widetilde{ x } }_{i} \right|}}


   The prices are scaled that they are unity in the base period or - if there
   is more than one base period - that the
   means of the prices over the base periods are unity.
   The argument \code{base} can be either \cr
   (a) a single number: the row number of the base prices, \cr
   (b) a vector indicating several observations: The means of these
   observations are used as base prices, \cr
   (c) a logical vector with the same length as the \code{data}: The
   means of the observations indicated as 'TRUE' are used as base prices, or
   (d) \code{NULL}: prices are not scaled.
}

\seealso{\code{\link{snqProfitEst}}.}

\author{Arne Henningsen}

\examples{
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput   / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   priceNames <- c( "pOutput", "pVarInput", "pLabor" )
   quantNames <- c( "qOutput", "qVarInput", "qLabor" )
   snqProfitWeights( priceNames, quantNames, germanFarms )
}

\keyword{models}

