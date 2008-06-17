.quadFuncCoefNames <- function( nExog ) {
   result <- paste( "alpha", c( 0:nExog ), sep = "_" )
   for( i in 1:nExog ) {
      for( j in i:nExog ) {
         result <- c( result, paste( "beta", i, j, sep = "_" ) )
      }
   }
   return( result )
}
