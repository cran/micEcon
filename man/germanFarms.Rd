\name{germanFarms}
\docType{data}
\alias{germanFarms}

\title{Output and Inputs of Farms in West-Germany}

\description{
   The \code{germanFarms} data frame contains annual data of an average
   full-time farm in West-Germany. Additionally, the price indices
   for agricultural output and agricultural variable input are
   included. 20 book-keeping years are included - starting in
   1975/76 and ending in 1994/95.
}

\usage{data(germanFarms)}

\format{
   This data frame contains the following columns:
   \describe{
      \item{year}{the book-keeping year.}
      \item{vOutput}{the value of Outputs (in current Deutschmark).}
      \item{pOutput}{price index of agricultural outputs (1980/81 = 100).}
      \item{vVarInput}{the value of variable inputs (in current Deutschmark).}
      \item{pVarInput}{price index of variable agricultural inputs (1980/81 = 100).}
      \item{qLabor}{the number of full-time worker equivalents.}
      \item{pLabor}{costs of an agricultural worker (Deutschmarks per year).}
      \item{land}{land used for agricultural production (in ha).}
   }
}

\source{
   Bundesministerium f\"{u}r Ern\"{a}hrung, Landwirtschaft und Forsten
   (Federal Department for Food, Agriculture and Forests),
   \emph{Agrarbericht der Bundesregierung
   (Agricultural Report of the Federal Government)},
   Jahrg\"{a}nge 1977-1996
   (years 1977-1996).
}

\keyword{datasets}
