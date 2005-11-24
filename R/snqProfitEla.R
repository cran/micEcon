## ===== calculation of elasticities from beta matrix ===
snqProfitEla <- function( beta, prices, quant, weights ) {
   if( !is.matrix( beta ) ) {
      stop( "argument 'beta' must be a matrix" )
   }
   if( nrow( beta ) != ncol( beta ) ) {
      stop( "argument 'beta' must be a quadratic matrix" )
   }
   if( length( prices ) != length( quant ) ) {
      stop( "arguments 'prices' and 'quant' must have the same length" )
   }
   if( length( prices ) != length( weights ) ) {
      stop( "arguments 'prices' and 'weights' must have the same length" )
   }
   if( nrow( beta ) != length( prices ) ) {
      stop( "arguments 'prices' must have as many elements as",
         " argument 'beta' has rows" )
   }
   nNetput  <- ncol( beta )
   prices   <- unlist( prices )
   quant    <- unlist( quant )
   hessian  <- snqProfitHessian( beta, prices, weights )
   ela      <- hessian * array( 1, c( nNetput ) ) %*% t( prices ) /
                  quant %*% t( array( 1, c( nNetput ) ) )
   if( !is.null( names( quant ) ) ) {
      rownames( ela ) <- names( quant )
   } else {
      rownames( ela ) <- paste( "q", 1:nNetput, sep = "" )
   }
   if( !is.null( names( prices ) ) ) {
      colnames( ela ) <- names( prices )
   } else {
      colnames( ela ) <- paste( "p", 1:nNetput, sep = "" )
   }
   return( ela )
}
