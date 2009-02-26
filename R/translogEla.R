translogEla <- function( xNames, data, coef, coefCov = NULL,
   quadHalf = TRUE, dataLogged = FALSE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef > length( coef ) ) {
      stop( "a translog function with ", nExog, " exogenous variables",
         " must have at least ", nCoef, " coefficients" )
   }

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- .micEconLogData( data = data, 
         varNames = xNames )
   }

   result <- quadFuncDeriv( xNames = xNames, data = logData, coef = coef, 
      coefCov = coefCov, quadHalf = quadHalf )

   names( result )[ 1 ] <- "ela"
   class( result )  <- "translogEla"

   return( result )
}
