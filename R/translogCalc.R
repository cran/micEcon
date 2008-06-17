translogCalc <- function( xNames, data, coef, quadHalf = TRUE,
   dataLogged = FALSE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef != length( coef ) ) {
      stop( "a translog function with ", nExog, " exogenous variables",
         " must have exactly ", nCoef, " coefficients" )
   }

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- data.frame( no = c( 1:nrow( data ) ) )
      for( i in seq( along = xNames ) ) {
         logData[[ xNames[ i ] ]] <- log( data[[ xNames[ i ] ]] )
      }
   }

   result <- quadFuncCalc( xNames, logData, coef, quadHalf = quadHalf )

   if( !dataLogged ) {
      result <- exp( result )
   }

   return( result )
}
