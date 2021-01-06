\name{appleProdFr86}
\docType{data}
\alias{appleProdFr86}

\title{Data on French Apple Producers in 1986}

\description{
   The \code{appleProdFr86} data frame includes 
   cross-sectional production data of 140 French apple producers 
   from the year 1986. 
   These data have been extracted from a panel data set
   that was used in Ivaldi et al. (1996).

}

\usage{data(appleProdFr86)}

\format{
   This data frame contains the following columns:
   \describe{
      \item{vCap}{costs of capital (including land).}
      \item{vLab}{costs of labour 
         (including remuneration of unpaid family labour).}
      \item{vMat}{costs of intermediate materials 
         (e.g. seedlings, fertilizer, pesticides, fuel).}
      \item{qApples}{quantity index of produced apples.}
      \item{qOtherOut}{quantity index of all other outputs.}
      \item{qOut}{quantity index of all outputs
         (not in the original data set, 
         calculated as 580,000 \eqn{\cdot}{*} 
         (\code{qApples} + \code{qOtherOut})).}
      \item{pCap}{price index of capital goods}
      \item{pLab}{price index of labour.}
      \item{pMat}{price index of materials.}
      \item{pOut}{price index of the aggregate output 
         (not in the original data set, artificially generated).}
      \item{adv}{dummy variable indicating the use of an advisory service
         (not in the original data set, artificially generated).}
   }
}

\source{
   This cross-sectional data set has been extracted 
   from a panel data set
   that is available in the data archive 
   of the Journal of Applied Econometrics: 
   \url{http://qed.econ.queensu.ca/jae/1996-v11.6/ivaldi-ladoux-ossard-simioni/}. 
}

\references{
   Ivaldi, M., N. Ladoux, H. Ossard, and M. Simioni (1996)
   Comparing Fourier and Translog Specifications of Multiproduct 
   Technology: Evidence from an Incomplete Panel of French Farmers.
   \emph{Journal of Applied Econometrics}, 11(6), p. 649-667.
}

\keyword{datasets}
