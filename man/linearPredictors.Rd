\name{linearPredictors}
\alias{linearPredictors}
\alias{linearPredictors.probit}

\title{Calculates linear predictors for different models}

\description{
  Calculates the (unobservable) linear predictors for probability models.
}

\usage{
linearPredictors(x, ...)
}

\arguments{
  \item{x}{model of an appropriate class}
  \item{\dots}{other arguments depending on the method}
}

\details{
  It is a generic function with a method for 'probit'.
}

\value{
  A matrix with nrow equal to the number of observations and one column:
  the linear predictors for observations
}

\author{Ott Toomet \email{otoomet@ut.ee}}

\examples{
data(Mroz87)
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
a <- probit(lfp ~ kids + educ + hushrs + huseduc + huswage + mtr +
 motheduc, data=Mroz87)
b <- linearPredictors(a)
cor(Mroz87$lfp, b) # should be positive and highly significant
}
\keyword{models}
\keyword{methods}
