aidsTestConsist <- function( priceNames, shareNames, totExpName, data = NULL,
   coef = NULL, alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) ) {

   if( length( priceNames ) != length( shareNames ) ) {
      stop( "arguments 'priceNames' and 'shareNames' must have the same length" )
   }

   result <- list()
   nGoods <- length( priceNames )
   nObs <- nrow( data )

   xt <- data[[ totExpName ]]
   priceMat <- array( NA, c( nObs, nGoods ) )
   shareMat <- array( NA, c( nObs, nGoods ) )
   for( i in 1: nGoods ) {
      priceMat[ , i ] <- data[[ priceNames[ i ] ]]
      shareMat[ , i ] <- data[[ shareNames[ i ] ]]
   }
   fitted <- aidsCalc( priceNames, totExpName, data, alpha0 = alpha0, coef = coef )

   # testing for monotonicity
   mono <- array( TRUE, c( nObs ) )
   cMatrices <- list()    # testing for concavity
   conc <- array( TRUE, c( nObs ) )

   lnp <- aidsPx( "TL", priceNames, data = data,
      alpha0 = alpha0, coef = coef )

   for( t in 1:nObs ) {
      mono[ t ] <- ( min( fitted$shares[ t, ] ) >= 0 )
      cMatrices[[ t ]] <- coef$gamma + ( coef$beta %*% t( coef$beta ) ) *
         ( log( xt[ t ] ) - lnp[ t ] ) -
         diag( shareMat[ t, ] ) + shareMat[ t, ] %*% t( shareMat[ t, ] )

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
