\name{rSquared}
\alias{rSquared}
\title{Calculate R squared value}

\description{
   Calculate R squared value.
}

\usage{rSquared( y, resid )}

\arguments{
   \item{ y }{vector of endogenous variables}
   \item{ resid }{vector of residuals}
}

\author{Arne Henningsen}

\examples{
   data( Blanciforti86 )
   reg <- lm( wFood1 ~ pFood1 + xFood, Blanciforti86 )
   rSquared( Blanciforti86$wFood1, reg$residuals )
   summary( reg )$r.squared  # returns the same value
}

\keyword{univar}
\keyword{multivariate}
\keyword{array}

