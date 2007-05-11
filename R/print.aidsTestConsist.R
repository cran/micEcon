print.aidsTestConsist <- function( x, ... ) {
   cat( "\nTests for theoretical consistency of an estimated " )
   cat( "Almost Ideal Demand System (AIDS):\n" )
   cat( "Observation where monotonicity is fulfilled: " )
   cat( x$mPercent, "%\n" )
   cat( "Observation where concavity is fulfilled: " )
   cat( x$cPercent, "%\n" )
}
