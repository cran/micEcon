\name{semidefiniteness}
\alias{semidefiniteness}
\title{Test for negative and positive semidefiniteness}

\description{
   Test wether a symmetric matrix is negative semidefinite,
   positive semidefinite, both of them or none of them.
}

\usage{semidefiniteness( m, tol = .Machine$double.eps, method = "det" )}

\arguments{
   \item{m}{a quadratic matrix}
   \item{tol}{tolerance level (values between \code{-tol} and \code{tol}
      are considered to be zero).}
   \item{method}{method to test for semidefiniteness, either "det"
      (the textbook method: checking for the signs of the determinants
      of sub-matrices) or "eigen" (checking for the signs of the eigen values).}
}

\references{
   Chiang, A.C. (1984)
   \emph{Fundamental Methods of Mathematical Economics},
   3rd ed., McGraw-Hill.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
   # a positive semidefinite matrix
   semidefiniteness( matrix( 1, 3, 3 ))

   # a negative semidefinite matrix
   semidefiniteness( matrix(-1, 3, 3 ))

   # a matrix that is positive and negative semidefinite
   semidefiniteness( matrix( 0, 3, 3 ))

   # a matrix that is neither positive nor negative semidefinite
   semidefiniteness( matrix( 1:9, 3, 3 ))
}

\keyword{array}
