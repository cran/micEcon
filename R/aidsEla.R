aidsEla <- function( coef, W, P = NULL, formula = "AIDS",
   qNames = NULL, pNames = NULL ) {

   if( length( coef$alpha ) != length( coef$beta ) ) {
      stop( "Arguments 'alpha' and 'beta' must have the same length." )
   } else if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
      stop( "Argument 'gamma' must be a square matrix" )
   } else if( length( coef$alpha ) != nrow( coef$gamma ) ) {
      stop( paste( "Number of rows of argument 'gamma' must be equal",
         "to the length of argument 'alpha'" ) )
   } else if(  length( coef$alpha ) != length( W ) ) {
      stop( "Arguments 'alpha' and 'W' must have the same length." )
   } else if(  length( coef$alpha ) != length( P ) && !is.null( P ) ) {
      stop( "Arguments 'alpha' and 'P' must have the same length." )
   }

   if( formula %in% c( "AIDS" ) ) {
      if( is.null( P ) ) {
         stop( "The 'AIDS' formula requires argument 'P' (prices)." )
      }
   } else if( formula %in% c( "Ch", "EU" ) ) {
      if( !is.null( P ) ) {
         warning( "The 'Ch' and 'EU' formulas do not require argument 'P' (prices)." )
      }
   }

   nGoods <- length( coef$alpha )

   ela <- list()
   ela$formula <- formula

   W <- array(W)

   if( formula == "AIDS" ) {
      P <- array(P)
      ela$exp <- array( 1, c( nGoods ) ) + coef$beta/W
      ela$hicks <- -diag( 1, nGoods, nGoods ) +
         array( 1, c( nGoods )) %*% t( W ) +
         coef$gamma / ( W %*% t( array( 1, c( nGoods )))) -
         coef$beta %*% t( array( 1, c( nGoods ))) *
         ( array( 1, c( nGoods )) %*% t( coef$alpha ) -
         array( 1, c( nGoods )) %*% t( W )+
         array( 1, c( nGoods )) %*% t( coef$gamma %*% log( P ))) /
         ( W %*% t( array( 1, c( nGoods ))))
      ela$marshall <- ela$hicks - ( ela$exp %*% t( array( 1, c( nGoods )))) *
         ( array( 1, c( nGoods )) %*% t( W ))
   } else if(formula=="Ch") {
      ela$exp <- array( 1, c( nGoods ) ) + coef$beta / W
      ela$hicks <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( W %*% t( array( 1, c( nGoods ) ) ) ) +
         array( 1, c( nGoods )) %*% t( W )
      ela$marshall <- ela$hicks - ( ela$exp %*% t( array( 1, c( nGoods ) ) ) ) *
         ( array( 1, c( nGoods )) %*% t(W))
   } else if( formula == "EU" ) {
      ela$exp <- array( 1, c( nGoods ) ) + coef$beta / W
      ela$marshall <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( W %*% t( array( 1, c( nGoods ) ) ) )
      ela$hicks <- ela$marshall + ( ela$exp %*% t( array( 1, c( nGoods ) ) ) ) *
         ( array( 1, c( nGoods )) %*% t(W))
   } else if( formula == "B1 not implemente" ) {
      ela$exp <- array( 1, c( nGoods ) ) + coef$beta / W
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / W[ i ] -
               coef$beta[ i ] * W[ j ] / W[ i ]
            for( k in 1:nGoods ) {
               ela$marshall[ i, j ] <- ela$marshall[ i, j ] - coef$beta[ i ] *
                  W[ k ] * log( P[ k ] ) * ( 0 + ( k == j ) ) / W[ i ]
            }
         }
      }
      ela$hicks <- ela$marshall + ( ela$exp %*% t( array( 1, c( nGoods ) ) ) ) *
         ( array( 1, c( nGoods )) %*% t(W))
   } else {
      stop( paste( "Formula '", as.character( formula ), "' is not supported.",
         sep = "" ) )
   }
   if( !is.null( qNames ) ) {
      rownames( ela$hicks ) <- qNames
      rownames( ela$marshall ) <- qNames
      names( ela$exp ) <- qNames
   }
   if( !is.null( pNames ) ) {
      colnames( ela$hicks ) <- pNames
      colnames( ela$marshall ) <- pNames
   }
   return( ela )
}
