\name{colMedians}
\alias{colMedians}
\alias{rowMedians}
\title{Medians of Columns and Rows}

\description{
   Compute the sample medians of the columns or rows of a data.frame or matrix.
}

\usage{
colMedians( x, na.rm = FALSE )

rowMedians( x, na.rm = FALSE )
}

\arguments{
   \item{x}{a data.frame or matrix.}
   \item{na.rm}{a logical value indicating whether \code{NA} values should be
      stripped before the computation proceeds.}
}

\value{
   A vector of the medians of each column or row of \code{x}.
}

\seealso{\code{\link{median}},\code{\link{colMeans}}.}

\author{Arne Henningsen}

\examples{
   data( "germanFarms" )
   colMedians( germanFarms[ , -1 ] )
}

\keyword{array}
