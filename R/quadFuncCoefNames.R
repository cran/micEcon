.quadFuncCoefNames <- function( nExog, nShifter = 0 ) {
   result <- paste( "alpha", c( 0:nExog ), sep = "_" )
   if( nExog > 0 ) {
      for( i in 1:nExog ) {
         for( j in i:nExog ) {
            result <- c( result, paste( "beta", i, j, sep = "_" ) )
         }
      }
   }
   if( nShifter > 0 ){
      result <- c( result, paste( "delta", 1:nShifter, sep = "_" ) )
   }
   return( result )
}
