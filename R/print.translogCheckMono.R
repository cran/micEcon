print.translogCheckMono <- function( x, ... ) {

   cat( "This translog function is" )
   if( x$strict ) {
      cat( " strictly" )
   }
   cat( " monotonically" )
   if( x$increasing ) {
      cat( " increasing" )
   } else {
      cat( " decreasing" )
   }
   cat( " at", sum( x$obs, na.rm = TRUE ) )
   cat( " out of", sum( !is.na( x$obs ) ), "observations (" )
   cat( round( 100 * sum( x$obs, na.rm = TRUE ) / sum( !is.na( x$obs ) ),
      digits = 1 ) )
   cat( "%)\n" )
   invisible( x )
}
