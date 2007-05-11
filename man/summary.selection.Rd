\name{summary.selection}
\alias{summary.selection}
\alias{print.summary.selection}

\title{Summarizing Selection Estimations}

\description{
   Print or return a summary of a selection estimation.
}

\usage{
   \method{summary}{selection}(object, ...)
   \method{print}{summary.selection}(x,
      digits = max(3, getOption("digits") - 3),
      part = "full", ...)
}

\arguments{
   \item{object}{an object of class '\code{selection}'.}
   \item{x}{an object of class '\code{summary.selection}'.}
   \item{part}{character string: which parts of the summary to print:
      "full" for all the estimated parameters
      (probit selection, outcome estimates, correlation and residual variance),
      or "outcome" for the outcome results only.}
   \item{digits}{numeric, (suggested) number of significant digits.}
   \item{...}{currently not used.}
}

\details{
   The variance-covariance matrix of the two-step estimator is currently
   implemented only for tobit-2 (sample selection) models,
   but not for the tobit-5 (switching regression) model.
}

\value{
   Summary methods return an object of class \code{summary.selection}.
   Print methods return the argument invisibly.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}, Ott
  Toomet \email{otoomet@ut.ee}}

\seealso{\code{\link{selection}}}

\examples{
## Wooldridge( 2003 ): example 17.5, page 590
data( Mroz87 )
wooldridge <- selection( lfp ~ nwifeinc + educ + exper + I( exper^2 ) +
   age + kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ),
   data = Mroz87, method = "2step" )

# summary of the 1st step probit estimation (Example 17.1, p. 562f)
# and the 2nd step OLS regression (example 17.5, page 590)
summary( wooldridge )

# summary of the outcome equation only
print(summary(wooldridge), part="outcome")
}

\keyword{models}



