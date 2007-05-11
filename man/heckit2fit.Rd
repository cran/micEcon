\name{heckit2fit}
\alias{heckit2fit}
\alias{heckit5fit}

\title{2-step Heckman (heckit) estimation}

\description{
   These functions do the actual fitting of tobit-2
   (sample selection) and tobit-5 (switching regression)
   models by the 2-step Heckman (heckit) estimation.
   They are called by \code{\link{selection}} or
   \code{\link{heckit}} and
   they are intended for \code{micEcon} internal use.
}

\usage{
heckit2fit( selection, outcome, data, inst = NULL, print.level = 0, ... )

heckit5fit( selection, outcome1, outcome2, data = sys.frame(sys.parent()),
   ys = FALSE, yo = FALSE, xs = FALSE, xo = FALSE, mfs = FALSE,
   mfo = FALSE, print.level = 0, ... )
}


\arguments{
  \item{selection}{formula for the probit estimation (1st step) (see details).}
  \item{outcome}{formula to be estimated (2nd step).}
  \item{outcome1}{formula, the first outcome equation.}
  \item{outcome2}{formula, the second outcome equation.}
  \item{data}{a data frame containing the data.}
  \item{inst}{an optional one-sided formula specifying instrumental variables
    for a 2SLS/IV estimation on the 2nd step.}
  \item{ys, yo, xs, xo, mfs, mfo}{logicals.  If true, the response (\code{y}),
    model matrix (\code{x}) or the model frame (\code{mf})
    of the selection (\code{s}) or outcome
    (\code{o}) equation(s) are returned.}
  \item{print.level}{numeric, values greater than 0 will produce
    increasingly more debugging information.}
  \item{...}{currently not used.}
}

\value{
   see \code{\link{selection}}.
}

\references{
   see \code{\link{selection}}.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de},
   Ott Toomet \email{otoomet@ut.ee}}

\seealso{\code{\link{selection}}, \code{\link{heckit}}}

\keyword{internal}
\keyword{models}
\keyword{regression}
