aidsTestConsist <- function( pNames, wNames, xtName, data = NULL,
   coef = NULL, alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) ) {

   if( length( pNames ) != length( wNames ) ) {
      stop( "arguments 'pNames' and 'wNames' must have the same length" )
   }

   result <- list()
   nGoods <- length( pNames )
   nObs <- nrow( data )

   xt <- data[[ xtName ]]
   prices <- array( NA, c( nObs, nGoods ) )
   shares <- array( NA, c( nObs, nGoods ) )
   for( i in 1: nGoods ) {
      prices[ , i ] <- data[[ pNames[ i ] ]]
      shares[ , i ] <- data[[ wNames[ i ] ]]
   }
   fitted <- aidsCalc( pNames, xtName, data, alpha0 = alpha0, coef = coef )

   # testing for monotonicity
   mono <- array( TRUE, c( nObs ) )
   cMatrices <- list()    # testing for concavity
   conc <- array( TRUE, c( nObs ) )

   lnp <- aidsPx( "TL", pNames, data = data,
      alpha0 = alpha0, coef = coef )

   for( t in 1:nObs ) {
      mono[ t ] <- ( min( fitted$shares[ t, ] ) >= 0 )
      cMatrices[[ t ]] <- coef$gamma + ( coef$beta %*% t( coef$beta ) ) *
         ( log( xt[ t ] ) - lnp[ t ] ) -
         diag( shares[ t, ] ) + shares[ t, ] %*% t( shares[ t, ] )

  #    for( i in 1:nGoods ) {
  #       conc[ t ] <- ( conc[ t ] & cMatrices[[ t ]][ i, i ] <= 0 )
  #    }
  #    for( i in 2:( nGoods - 1 ) ) {
  #       conc[ t ] <- ( conc[ t ] &
  #      ( det( cMatrices[[ t ]][ 1:i, 1:i ] ) * (-1)^i >= 0 ) )
  #    }
      conc[ t ] <- semidefiniteness( cMatrices[[ t ]][ 1:( nGoods - 1),
         1:( nGoods - 1) ] )$negative
   }
   result$mPercent <- 100 * sum( mono ) / nObs
   result$monotony <- mono
   result$cPercent <- 100 * sum( conc ) / nObs
   result$concavity <- conc
   result$cMatrices <- cMatrices
   class( result ) <- "aidsTestConsist"
   return( result )
}
