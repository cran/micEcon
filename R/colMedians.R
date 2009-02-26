colMedians <- function( x, na.rm = FALSE ) {

   if( is.data.frame( x ) ) {
      x <- as.matrix( x )
   }
   if( !is.matrix( x ) ) {
      stop( "argument 'x' must be a data.frame or a matrix" )
   }
   if( !is.numeric( x ) ) {
      stop( "argument 'x' must be numeric" )
   }

   result <- rep( NA, ncol( x ) )
   names( result ) <- colnames( x )

   for( i in 1:ncol( x ) ) {
      result[ i ] <- median( x[ , i ], na.rm = na.rm )
   }

   return( result )
}

rowMedians <- function( x, na.rm = FALSE ) {
   return( colMedians( t( x ), na.rm = na.rm ) )
}
