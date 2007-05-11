\name{predict.snqProfitEst}
\alias{predict.snqProfitEst}
\alias{predict.snqProfitImposeConvexity}

\title{Predictions from an SNQ profit function}

\description{
   Returns the predicted values, their standard errors and the
   confidence limits of prediction for an
   Symmetric Normalized Quadratic (SNQ) profit function.
}

\usage{
\method{predict}{snqProfitEst}( object, newdata = object$data,
   se.fit = FALSE, se.pred = FALSE, interval = "none", level = 0.95, ... )

\method{predict}{snqProfitImposeConvexity}( object, newdata = object$data,
   se.fit = FALSE, se.pred = FALSE, interval = "none", level = 0.95, ... )
}

\arguments{
  \item{object}{an object of type \code{snqProfitEst} or
      \code{snqProfitImposeConvexity}.}
  \item{newdata}{data frame in which to predict.}
  \item{se.fit}{logical. Return the standard error of the fitted values?}
  \item{se.pred}{logical. Return the standard error of prediction?}
  \item{interval}{Type of interval calculation ("none", "confidence"
      or "prediction").}
  \item{level}{confidence level.}
  \item{...}{currently not used.}
}

\details{
   The variance of the fitted values
   (used to calculate the standard errors of the fitted values
   and the "confidence interval") is calculated by
   \eqn{Var[E[y^0]-\hat{y}^0]=x^0 \; Var[b] \; {x^0}'}\cr
   The variances of the predicted values
   (used to calculate the standard errors of the predicted values
   and the "prediction intervals") is calculated by
   \eqn{Var[y^0-\hat{y}^0]=\hat{\sigma}^2+x^0 \; Var[b] \; {x^0}'}
}

\value{
   \code{predict.snqProfitEst} and \code{predict.snqProfitImposeConvexity} return
   a dataframe that
   contains the predicted profit and for each netput the predicted quantities
   (e.g. "quant1" ) and if requested
   the standard errors of the fitted values (e.g. "quant1.se.fit"),
   the standard errors of the prediction (e.g. "quant1.se.pred"),
   and the lower (e.g. "quant1.lwr") and upper (e.g. "quant1.upr")
   limits of the confidence or prediction interval(s).
}

\references{
   Diewert, W.E. and T.J. Wales (1987)
   Flexible functional forms and global curvature conditions.
   \emph{Econometrica}, 55, p. 43-68.

   Diewert, W.E. and T.J. Wales (1992)
   Quadratic Spline Models for Producer's Supply and Demand Functions.
   \emph{International Economic Review}, 33, p. 705-722.

   Greene, W. H. (2003)
   \emph{Econometric Analysis, Fifth Edition}, Macmillan.

   Gujarati, D. N. (1995)
   \emph{Basic Econometrics, Third Edition}, McGraw-Hill.

   Kmenta, J. (1997)
   \emph{Elements of Econometrics, Second Edition}, University of
   Michigan Publishing.

   Kohli, U.R. (1993)
   A symmetric normalized quadratic GNP function and the US demand
   for imports and supply of exports.
   \emph{International Economic Review}, 34, p. 243-255.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}


\seealso{
   \code{\link{snqProfitEst}}, \code{\link{snqProfitCalc}}
   and \code{\link{predict}}
}

\examples{
   data( germanFarms )
   germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
   germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
   germanFarms$qLabor    <- -germanFarms$qLabor
   germanFarms$time      <- c( 0:19 )
   priceNames <- c( "pOutput", "pVarInput", "pLabor" )
   quantNames <- c( "qOutput", "qVarInput", "qLabor" )
   estResult <- snqProfitEst( priceNames, quantNames, c("land","time"), data=germanFarms )
   predict( estResult )
   predict( estResult, se.fit = TRUE, se.pred = TRUE, interval = "confidence" )
}

\keyword{models}



