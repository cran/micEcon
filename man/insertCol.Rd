\name{insertCol}
\alias{insertCol}
\title{Insert Column into a Matrix}

\description{
   Insert Column into a Matrix.
}

\usage{insertCol( m, c, v = NA )}

\arguments{
   \item{ m }{Matrix.}
   \item{ c }{column number where the new column should be inserted.}
   \item{ v }{Optional values of the new column.}
}

\value{
   a matrix with one more column than the provided matrix \code{m}.
}

\seealso{\code{\link{insertRow}}.}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   m <- matrix( 1:4, 2 )
   insertCol( m, 2, 5:6 )
}

\keyword{array}
