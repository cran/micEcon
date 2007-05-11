.aidsEstMethod <- function( method, priceIndex ) {

   if( substr( method, 1, 2 ) == "LA" ) {
      result <- "Linear Approximation (LA) with"
      if( priceIndex == "S" ) {
         result <- paste( result, "Stone Index (S)\n" )
      } else if( priceIndex == "SL" ) {
         result <- paste( result, "lagged Stone Index (SL)\n" )
      } else if( priceIndex == "P" ) {
         result <- paste( result, "Paasche Index (P)\n" )
      } else if( priceIndex == "L" ) {
         result <- paste( result, "Laspeyres Index (L)\n" )
      } else if( priceIndex == "T" ) {
         result <- paste( result, "Tornqvist Index (T)\n" )
      } else {
         result <- paste( result, "unknown price index\n" )
      }
   } else if( substr( method, 1, 2 ) %in% c( "MK", "IL" ) ) {
      result <- "'Iterated Linear Least Squares Estimator' (IL)\n(starting with"
      if( priceIndex == "S" ) {
         result <- paste( result, "Stone Index, S)\n" )
      } else if( priceIndex == "SL" ) {
         result <- paste( result, "lagged Stone Index, SL)\n" )
      } else if( priceIndex == "P" ) {
         result <- paste( result, "Paasche Index, P)\n" )
      } else if( priceIndex == "L" ) {
         result <- paste( result, "Laspeyres Index, L)\n" )
      } else if( priceIndex == "T" ) {
         result <- paste( result, "Tornqvist Index, T)\n" )
      } else {
         result <- paste( result, "unknown price index)\n" )
      }
   } else {
      result <- "unknown method"
   }

   return( result )
}
