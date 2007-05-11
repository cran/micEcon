aidsCalc <- function( priceNames, totExpName, data = NULL, px = "TL", lnp = NULL,
   coef = NULL, alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) ) {

   if( px != "TL" && is.null( lnp ) ) {
      stop( "at the moment only the translog (TL) price index works",
         " if argument 'lnp' is not specified" )
   }
   nGoods <- length( priceNames )

   if( is.null( lnp ) ) {
      lnp <- aidsPx( px, priceNames, data = data,
         alpha0 = alpha0, coef = coef )
   }
   shareData <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( shareData ) <- paste( "w", as.character( 1:nGoods ), sep = "" )
   rownames( shareData ) <- rownames( data )
   quant <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( quant ) <- paste( "q", as.character( 1:nGoods ), sep = "" )
   rownames( quant ) <- rownames( data )
   for( i in 1:nGoods ) {
      shareData[ , i ] <- coef$alpha[ i ] + coef$beta[ i ] *
         ( log( data[[ totExpName ]] ) - lnp )
      for( j in 1:nGoods ) {
         shareData[ , i ] <- shareData[ , i ] + coef$gamma[ i, j ] *
            log( data[[ priceNames[ j ] ]] )
      }
      quant[ , i ] <- shareData[ , i ] * data[[ totExpName ]] / data[[ priceNames[ i ] ]]
   }
   result <- list()
   result$shares <- shareData
   result$quant  <- quant
   return( result )
}
