\name{insertRow}
\alias{insertRow}
\title{Insert Row into a Matrix}

\description{
   Insert Row into a Matrix.
}

\usage{insertRow( m, r, v = NA )}

\arguments{
   \item{ m }{Matrix.}
   \item{ r }{Row number where the new Rown should be inserted.}
   \item{ v }{Optional values for the new Row.}
}

\value{
   a matrix with one more row than the provided matrix \code{m}.
}

\seealso{\code{\link{insertCol}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   m <- matrix( 1:4, 2 )
   insertRow( m, 2, 5:6 )
}

\keyword{array}
