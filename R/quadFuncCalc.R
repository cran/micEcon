quadFuncCalc <- function( xNames, data, allCoef, quadHalf = TRUE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef != length( allCoef ) ) {
      stop( paste( "A quadratic function with", nExog, "exogenous variables",
         "must have exactly", nCoef, "coefficients." ) )
   }

   alpha0 <- allCoef[ 1 ]
   alpha  <- allCoef[ 2:( nExog + 1 ) ]
   beta   <- vecli2m( allCoef[ ( nExog + 2 ):nCoef ] )

   result <- rep( alpha0, nrow( data ) )
   for( i in 1:nExog ) {
      result <- result + alpha[ i ] * data[[ xNames[ i ] ]]
      for( j in 1:nExog ) {
         result <- result + ifelse( quadHalf, 0.5, 1 ) * beta[ i, j ] *
            data[[ xNames[ i ] ]]  * data[[ xNames[ j ] ]]
      }
   }

   names( result ) <- rownames( data )
   return( result )
}
