aidsEla <- function( coef, W, P = NULL, formula = "AIDS",
   qNames = NULL, pNames = NULL, coefVcov = NULL, df = NULL ) {

   nGoods <- length( coef$alpha )

   if( length( coef$alpha ) != length( coef$beta ) ) {
      stop( "arguments 'alpha' and 'beta' must have the same length" )
   } else if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
      stop( "argument 'gamma' must be a square matrix" )
   } else if( length( coef$alpha ) != nrow( coef$gamma ) ) {
      stop( "number of rows of argument 'gamma' must be equal",
         " to the length of argument 'alpha'" )
   } else if(  length( coef$alpha ) != length( W ) ) {
      stop( "arguments 'alpha' and 'W' must have the same length" )
   } else if(  length( coef$alpha ) != length( P ) && !is.null( P ) ) {
      stop( "arguments 'alpha' and 'P' must have the same length" )
   }
   if( is.null( qNames ) ) {
      qNames <- .aidsQuantNames( W, coef, nGoods )
   } else {
      if( length( qNames ) != nGoods ) {
         stop( "argument 'qNames' must have ", nGoods, " elements" )
      }
   }
   if( is.null( pNames ) ) {
      pNames <- .aidsPriceNames( P, coef, nGoods )
   } else {
      if( length( pNames ) != nGoods ) {
         stop( "argument 'pNames' must have ", nGoods, " elements" )
      }
   }


   if( formula %in% c( "AIDS" ) ) {
      if( is.null( P ) ) {
         stop( "the 'AIDS' formula requires argument 'P' (prices)" )
      }
   } else if( formula %in% c( "Ch", "EU" ) ) {
      if( !is.null( P ) ) {
         warning( "the 'Ch' and 'EU' formulas do not require argument 'P' (prices)" )
      }
   }

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
   } else if( formula == "B1 not implemented" ) {
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
      stop( "formula '", as.character( formula ), "' is not supported" )
   }
   names( ela$exp )         <- qNames
   rownames( ela$hicks )    <- qNames
   colnames( ela$hicks )    <- pNames
   rownames( ela$marshall ) <- qNames
   colnames( ela$marshall ) <- pNames
   if( !is.null( coefVcov ) && formula %in% c( "AIDS" ) ) {
      jacobian <- aidsElaJacobian( coef = coef, share = W, price = P,
         formula = formula, qNames = qNames, pNames = pNames )
      ela$allVcov      <- jacobian$all      %*% coefVcov %*% t( jacobian$all )
      ela$expVcov      <- jacobian$exp      %*% coefVcov %*% t( jacobian$exp )
      ela$hicksVcov    <- jacobian$hicks    %*% coefVcov %*% t( jacobian$hicks )
      ela$marshallVcov <- jacobian$marshall %*% coefVcov %*% t( jacobian$marshall )
      # standard errors
      ela$expStEr      <- diag( ela$expVcov )^0.5
      ela$hicksStEr    <- matrix( diag( ela$hicksVcov )^0.5,
         ncol = nGoods, byrow = TRUE )
      ela$marshallStEr <-  matrix( diag( ela$marshallVcov )^0.5,
         ncol = nGoods, byrow = TRUE )
      # dim names for standard errors
      names( ela$expStEr )         <- names( ela$exp )
      dimnames( ela$hicksStEr )    <- dimnames( ela$hicks )
      dimnames( ela$marshallStEr ) <- dimnames( ela$marshall )
      # t-values
      ela$expTval      <- ela$exp      / ela$expStEr
      ela$hicksTval    <- ela$hicks    / ela$hicksStEr
      ela$marshallTval <- ela$marshall / ela$marshallStEr
      if( !is.null( df ) ) {
         ela$expPval <- 2 * pt( abs( ela$expTval ), df,
            lower.tail = FALSE )
         ela$hicksPval <- 2 * pt( abs( ela$hicksTval ), df,
            lower.tail = FALSE )
         ela$marshallPval <- 2 * pt( abs( ela$marshallTval ), df,
            lower.tail = FALSE )
      }
   }
   class( ela ) <- "aidsEla"
   return( ela )
}
