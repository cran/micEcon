quadFuncDeriv <- function( xNames, data, allCoef, allCoefCov = NULL,
   quadHalf = TRUE ) {

   checkNames( c( xNames ), names( data ) )

   result <- list()

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( nCoef != length( allCoef ) ) {
      stop( paste( "A quadratic function with", nExog, "exogenous variables",
         "must have exactly", nCoef, "coefficients." ) )
   }

   alpha  <- allCoef[ 2:( nExog + 1 ) ]
   beta   <- vecli2m( allCoef[ ( nExog + 2 ):nCoef ] )

   ## derivatives
   deriv <- array( NA, c( nrow( data ), nExog ) )
   for( i in 1:nExog ) {
      deriv[ , i ] <- alpha[ i ]
      for( j in 1:nExog ) {
         deriv[ , i ] <- deriv[ , i ] + ifelse( quadHalf, 1, 2 ) *
            beta[ i, j ] * data[[ xNames[ j ] ]]
      }
   }
   colnames( deriv ) <- xNames
   result$deriv    <- as.data.frame( deriv )

   if( !is.null( allCoefCov ) ) {
      ## variances of the derivatives
      variance <- array( NA, c( nrow( data ), nExog ) )
      for(i in 1:nExog ) {
         variance[ , i ] <- allCoefCov[ i + 1, i + 1 ]   # variance of aplha(i)
         for( j in 1:nExog ) {
            variance[ , i ] <- variance[ , i ] +
               allCoefCov[ i + 1, 1 + nExog + veclipos( i, j, nExog ) ] *
               ifelse( quadHalf, 1, 2 ) * data[[ xNames[ j ] ]]
               # covariance alpha(i)-beta(i,_)
         }
         for( j in 1:nExog ) {
            for( k in 1:nExog ) {
               variance[ , i ] <- variance[ , i ] +
                  allCoefCov[ 1 + nExog + veclipos( i, j, nExog ),
                  1 + nExog + veclipos( i, k, nExog ) ] *
                  ifelse( quadHalf, 1, 4 ) *
                  data[[ xNames[ j ] ]] * data[[ xNames[ k ] ]]
                  # variances + covariance beta(i,_)-beta(i,_)
            }
         }
      }
      stdDev <- variance^0.5  # standard errors
      colnames( variance ) <- xNames
      colnames( stdDev )   <- xNames
      result$variance <- as.data.frame( variance )
      result$stdDev   <- as.data.frame( stdDev )
   }

   class( result ) <- "quadFuncDeriv"
   return( result )
}
