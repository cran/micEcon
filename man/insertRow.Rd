\name{insertRow}
\alias{insertRow}
\title{Insert Row into a Matrix}

\description{
   Insert a new row into a matrix.
}

\usage{insertRow( m, r, v = NA )}

\arguments{
   \item{ m }{matrix.}
   \item{ r }{row number where the new row should be inserted.}
   \item{ v }{optional values for the new row.}
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
