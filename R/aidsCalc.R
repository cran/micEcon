aidsCalc <- function( pNames, xtName, data = NULL, px = "TL", lnp = NULL,
   coef = NULL, alpha0 = ifelse( is.null( coef$alpha0 ), 0, coef$alpha0 ) ) {

   if( px != "TL" && is.null( lnp ) ) {
      stop( paste( "At the moment only the translog (TL) price index works",
         "if argument 'lnp' is not specified" ) )
   }
   nGoods <- length( pNames )

   if( is.null( lnp ) ) {
      lnp <- aidsPx( px, pNames, data = data,
         alpha0 = alpha0, coef = coef )
   }
   shares <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( shares ) <- paste( "w", as.character( 1:nGoods ), sep = "" )
   rownames( shares ) <- rownames( data )
   quant <- as.data.frame( matrix( 0, nrow = nrow( data ), ncol = nGoods ) )
   names( quant ) <- paste( "q", as.character( 1:nGoods ), sep = "" )
   rownames( quant ) <- rownames( data )
   for( i in 1:nGoods ) {
      shares[ , i ] <- coef$alpha[ i ] + coef$beta[ i ] *
         ( log( data[[ xtName ]] ) - lnp )
      for( j in 1:nGoods ) {
         shares[ , i ] <- shares[ , i ] + coef$gamma[ i, j ] *
            log( data[[ pNames[ j ] ]] )
      }
      quant[ , i ] <- shares[ , i ] * data[[ xtName ]] / data[[ pNames[ i ] ]]
   }
   result <- list()
   result$shares <- shares
   result$quant  <- quant
   return( result )
}
