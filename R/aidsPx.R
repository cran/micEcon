aidsPx <- function( px, priceNames, shareNames = NULL, data = NULL, base = 1,
   coef = NULL, alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) ) {

   nGoods <- length( priceNames )
   if( !is.null( shareNames ) ) {
      if( nGoods != length( shareNames ) && px != "TL" ) {
         stop( "'shareNames' must have as many elements as 'priceNames'" )
      }
   }
   nObs <- nrow(  data )
   lnp <- array( 0, c( nObs ))
   if(px=="S") {      # Stone index
      for( i in 1:nGoods ) {
         lnp <- lnp + data[[ shareNames[ i ] ]] * log( data[[ priceNames[ i ] ]] )
      }
   } else if(px=="SL") {     # Stone index with lagged shares
      lnp[ 1 ] <- NA
      for( i in 1:nGoods ) {
         lnp[ 2:nObs ] <- lnp[ 2:nObs ] +
            data[[ shareNames[ i ] ]][ 1:(nObs-1) ] *
            log( data[[ priceNames[ i ] ]][ 2:nObs ] )
      }
   } else if(px=="P") {      # log-Paasche index
      for( i in 1:nGoods) {
         lnp <- lnp + data[[ shareNames[ i ] ]] * log( data[[ priceNames[ i ] ]] /
            mean( data[[ priceNames[ i ] ]][ base ] ) )
      }
   } else if(px=="L") {      # log-Laspeyres index
      for( i in 1:nGoods) {
         lnp <- lnp + mean( data[[ shareNames[ i ] ]][ base ] ) *
            log( data[[ priceNames[ i ] ]] )
      }
   } else if(px=="T") {      # Tornqvist index
      for( i in 1:nGoods) {
         lnp <- lnp + c( 0.5 * ( data[[ shareNames[ i ] ]] +
            mean( data[[ shareNames[ i ] ]][ base ] ) *
            matrix( 1, nrow = nObs ) ) * log( data[[ priceNames[ i ] ]] /
            mean( data[[ priceNames[ i ] ]][ base ] ) ) )
      }
   } else if(px=="TL") {      # Translog index
      lnp <- array( alpha0, c( nObs ) )
      for( i in 1:nGoods ) {
         lnp <- lnp + coef$alpha[ i ] * log( data[[ priceNames[ i ] ]] )
         for( j in 1:nGoods ) {
            lnp <- lnp + 0.5 * coef$gamma[ i, j ] *
               log( data[[ priceNames[ i ] ]] ) *
               log( data[[ priceNames[ j ] ]] )
         }
      }
   } else {
      stop( "the argument 'px' (price index) must be either 'S',",
         " 'SL', 'P', 'L', 'T' or 'TL'" )
   }
   return( lnp )
}
