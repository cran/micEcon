print.aidsEst <- function( x, ... ) {
   cat( "\nDemand analysis with the Almost Ideal " )
   cat( "Demand System (AIDS)\n" )
   cat( "\nEstimation Method: " )
   if( substr( x$method, 1, 2 ) == "LA" ) {
      cat( "Linear Approximation (LA) with " )
      if( x$px == "S" ) {
         cat( "Stone Index (S)\n" )
      } else if( x$px == "SL" ) {
         cat( "lagged Stone Index (SL)\n" )
      } else if( x$px == "P" ) {
         cat( "Paasche Index (P)\n" )
      } else if( x$px == "L" ) {
         cat( "Laspeyres Index (L)\n" )
      } else if( x$px == "T" ) {
         cat( "Tornqvist Index (T)\n" )
      } else {
         cat( "unknown price index\n" )
      }
   } else if( substr( x$method, 1, 2 ) == "MK" ) {
      cat( "Michalek & Keyzer (MK) starting with " )
      if( x$px == "S" ) {
         cat( "Stone Index (S)\n" )
      } else if( x$px == "SL" ) {
         cat( "lagged Stone Index (SL)\n" )
      } else if( x$px == "P" ) {
         cat( "Paasche Index (P)\n" )
      } else if( x$px == "L" ) {
         cat( "Laspeyres Index (L)\n" )
      } else if( x$px == "T" ) {
         cat( "Tornqvist Index (T)\n" )
      } else {
         cat( "unknown price index\n" )
      }
   }
   cat( "\nEstimated Coefficients\n" )
   cat( "alpha\n" )
   print( x$coef$alpha )
   cat( "beta\n" )
   print( x$coef$beta )
   cat( "gamma\n" )
   print( x$coef$gamma )

   cat( "\nDemand Elasticities " )
   if( x$ela$formula == "Ch" ) {
      cat( "(Formula of Chalfant / Goddard)\n" )
   } else if( x$ela$formula == "AIDS" ) {
      cat( "(original AIDS formula)\n" )
   }
   cat( "Expenditure Elasticities\n" )
   print( x$ela$exp )
   cat( "\nMarshallian (uncompensated) Price Elasticities\n" )
   print( x$ela$marshall )
   cat( "\nHicksian (compensated) Price Elasticities\n" )
   print( x$ela$hicks )
}
