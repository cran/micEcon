print.aidsElas <- function( x, ... ) {

   cat( "\nDemand Elasticities " )
   if( x$method %in% c( "GA", "B1" ) ) {
      cat( "(formulas of Green and Alston / Buse (first set))\n" )
   } else if( x$method %in% c( "B2" ) ) {
      cat( "(formulas of Buse (second set))\n" )
   } else if( x$method %in% c( "Go", "Ch" ) ) {
      cat( "(formulas of Goddard / Chalfant)\n" )
   } else if( x$method == "EU" ) {
      cat( "(formula of Eales and Unnevehr)\n" )
   } else if( x$method == "AIDS" ) {
      cat( "(original AIDS formulas)\n" )
   } else {
      cat( "(unknown formula '", x$method, "')\n", sep = "" )
   }
   cat( "Expenditure Elasticities\n" )
   print( x$exp )
   cat( "\nMarshallian (uncompensated) Price Elasticities\n" )
   print( x$marshall )
   cat( "\nHicksian (compensated) Price Elasticities\n" )
   print( x$hicks )
   invisible( x )
}
