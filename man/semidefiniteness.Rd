\name{semidefiniteness}
\alias{semidefiniteness}
\title{Positive or Negative Semidefiniteness}

\description{
   Check whether a symmetric matrix is positive or negative semidefinite.
}

\usage{semidefiniteness( m, positive = TRUE, tol = .Machine$double.eps,
   method = "det" )}

\arguments{
   \item{m}{a quadratic matrix or a list containing quadratic matrices.}
   \item{positive}{logical. Check for positive (\code{TRUE}, default)
      or negative (\code{FALSE}) semidefiniteness.}
   \item{tol}{tolerance level (values between \code{-tol} and \code{tol}
      are considered to be zero).}
   \item{method}{method to test for semidefiniteness, either "det"
      (the textbook method: checking for the signs of the determinants
      of sub-matrices) or "eigen" (checking for the signs of the eigen values).}
}

\details{
   Please note that a matrix can be
   neither positive nor negative semidefinite
   or positive and negative semidefinite
   at the same time.
}

\value{
   \code{semidefiniteness} returns a locigal value
   or a logical vector (if argument \code{m} is a list)
   indicating whether the matrix (or each of the matrices)
   is positive/negative (depending on argument \code{positive})
   semidefinite.
}

\references{
   Chiang, A.C. (1984)
   \emph{Fundamental Methods of Mathematical Economics},
   3rd ed., McGraw-Hill.
}

\author{Arne Henningsen}

\examples{
   # a positive semidefinite matrix
   semidefiniteness( matrix( 1, 3, 3 ))

   # a negative semidefinite matrix
   semidefiniteness( matrix(-1, 3, 3 ), positive = FALSE )

   # a matrix that is positive and negative semidefinite
   semidefiniteness( matrix( 0, 3, 3 ))
   semidefiniteness( matrix( 0, 3, 3 ), positive = FALSE )

   # a matrix that is neither positive nor negative semidefinite
   semidefiniteness( matrix( 1:9, 3, 3 ))
   semidefiniteness( matrix( 1:9, 3, 3 ), positive = FALSE )
}

\keyword{array}
