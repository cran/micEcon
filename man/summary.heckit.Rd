\name{summary.heckit}
\alias{summary.heckit}
\title{Print output of heckit estimation}

\description{
   This function returns a summary of the 2-step-Heckmann / heckit estimation.
}

\usage{
   \method{summary}{heckit}( object, ... )
}

\arguments{
   \item{object}{an object of calss \code{heckit}.}
   \item{...}{further arguments passed to \code{\link{summary.glm}}
      and \code{\link{summary.lm}}.}
}

\value{
  \code{summary.heckit} returns an object of class \code{heckit}.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\seealso{\code{\link{heckit}}}

\examples{
data( Mroz87 )
wooldridge <- heckit( log( wage ) ~ educ + exper + I( exper^2 ),
   lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age + kids5 + kids618, Mroz87 )
summary( wooldridge )        # summary of the 1st step probit estimation
                             # (Example 17.1, p. 562f) and 2nd step OLS regression
}

\keyword{models}



