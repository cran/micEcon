checkNames <- function( testNames, allNames ) {
   inAllNames   <- testNames %in% allNames
   if( !all( inAllNames ) ) {
      stop( paste( "Object(s) '",
         paste( testNames[ !inAllNames ], collapse = "', '" ),
         "' not found.", sep = "" ) )
   }
}

