\name{logDataSet}
\alias{logDataSet}
\title{Creating a Data Set with the Logarithms of the Original Variables}
\description{
   This function creates a data set with the logarithms
   of the original variables.
}

\usage{
logDataSet( data, varNames, varNamesNum = NULL  )
}

\arguments{
   \item{data}{a data frame containing the data
      (possibly a panel data frame created with
      \code{\link[plm]{pdata.frame}}).}
   \item{varNames}{vector of character strings
      that indicates names of variables in the data frame.
      The logarithm of these variables
      are included in the returned data frame.}
   \item{varNamesNum}{optional vector of character strings
      that indicates names of further variables in tha data frame.
      In case of numeric variables,
      the logarithms of these variabes
      are included in the returned data frame.
      In case of factor of logical variables,
      these variabes are included in the returned data frame
      without any transformation.}
}

\author{Arne Henningsen}

\examples{
   data( "germanFarms" )
   datLog <- logDataSet( germanFarms, c( "vAnimal", "vOutput", "vVarInput" ) )
   summary( datLog )
}

\keyword{array}
