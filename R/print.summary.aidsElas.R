print.summary.aidsElas <- function( x, ... ) {

   if( is.na( sum( x$table[ , 2 ] ) ) ) {
      print( x$table )
   } else {
      printCoefmat( x$table, ... )
   }

   invisible( x )
}