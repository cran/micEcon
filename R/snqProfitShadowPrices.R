snqProfitShadowPrices <- function( pNames, fNames, estResult = NULL,
   data = estResult$estData, weights = estResult$weights,
   coef = estResult$coef, form = estResult$form ) {

   checkNames( c( pNames, fNames ), names( data ) )

   nNetput <- length( pNames )
   nFix    <- length( fNames )
   nObs    <- nrow( data )

   snqProfitTestCoef( nNetput, nFix, coef, form = form,
      coefNames = c( "delta", "gamma" ) )

   normPrice <- numeric( nObs )
   for( i in 1:nNetput ) {
      normPrice <- normPrice + data[[ pNames[ i ] ]] * weights[ i ]
   }

   shadowPrices <- array( 0, c( nObs, nFix ) )
   for( j in 1:nFix ) {
      for( i in 1:nNetput ) {
         shadowPrices[ , j ] <- shadowPrices[ , j ] + coef$delta[ i, j ] *
            data[[ pNames[ i ] ]]
      }
      if( form == 0 ) {
         for( k in 1:nFix ) {
            shadowPrices[ , j ] <- shadowPrices[ , j ] + normPrice *
               coef$gamma[ j, k ] * data[[ fNames[ k ] ]]
         }
      } else {
         for( i in 1:nNetput ) {
            for( k in 1:nFix ) {
               shadowPrices[ , j ] <- shadowPrices[ , j ] +
                  coef$gamma[ i, j, k ] * data[[ pNames[ i ] ]] *
                  data[[ fNames[ k ] ]]
            }
         }
      }
   }

   result <- as.data.frame( shadowPrices )
   names( result ) <- fNames

   return( result )
}
