translogDeriv <- function( xNames, data, allCoef, allCoefCov = NULL,
   quadHalf = TRUE, logValues = FALSE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef != length( allCoef ) ) {
      stop( paste( "A translog function with", nExog, "exogenous variables",
         "must have exactly", nCoef, "coefficients." ) )
   }

   result <- list()

   alpha  <- allCoef[ 2:( nExog + 1 ) ]
   beta   <- vecli2m( allCoef[ ( nExog + 2 ):nCoef ] )

   if( logValues ) {
      logData <- data
   } else {
      logData   <- data.frame( no = c( 1:nrow( data ) ) )
      for( i in seq( along = xNames ) ) {
         logData[[ xNames[ i ] ]] <- log( data[[ xNames[ i ] ]] )
      }
   }

   logyHat <- translogCalc( xNames, logData, allCoef, quadHalf = quadHalf,
      logValues = TRUE )

   deriv <- matrix( 0, nrow( data ), nExog )
   for( i in seq( along = xNames ) ) {
      deriv[ , i ] <- alpha[ i ]
      for( j in seq( along = xNames ) ) {
         deriv[ , i ] <- deriv[ , i ] + ifelse( quadHalf, 1, 2 ) *
            beta[ i, j ] * logData[[ xNames[ j ] ]]
      }
      deriv[ , i ] <- deriv[ , i ] * exp( logyHat ) / exp( logData[[ xNames[ i ] ]] )
   }

   colnames( deriv ) <- xNames
   result$deriv      <- as.data.frame( deriv )

   class( result ) <- "translogDeriv"
   return( result )
}
