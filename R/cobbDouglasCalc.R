cobbDouglasCalc <- function( xNames, data, coef, dataLogged = FALSE) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )

   if( nExog + 1 != length( coef ) ) {
      stop( "a Cobb-Douglas function with ", nExog, " exogenous variables",
         " must have exactly ", nExog + 1, " coefficients" )
   }

   coefNames <- paste( "a", c( 0:nExog ), sep = "_" )
   if( is.null( names( coef ) ) ) {
      names( coef ) <- coefNames
   } else {
      coefMissing <- !( coefNames %in% names( coef ) )
      if( any( coefMissing ) ) {
         stop( "coefficient(s) ",
            paste( coefNames[ coefMissing ], collapse = ", " ),
            " are missing" )
      }
      rm( coefMissing )
   }
   rm( coefNames )

   if( dataLogged ) {
      logData <- data
   } else {
      logData <- .micEconLogData( data = data, varNames = xNames )
   }

   result <- rep( coef[ "a_0" ], nrow( data ) )
   for( i in seq( along = xNames ) ) {
      result <- result +
         coef[ paste( "a", i, sep = "_" ) ] * logData[[ xNames[ i ] ]]
   }

   if( !dataLogged ) {
      result <- exp( result )
   }

   names( result ) <- rownames( data )

   return( result )
}
