aidsJacobian <- function( allCoef, pNames, xtName, data = NULL,
      omitLast = TRUE, alpha0 = 0 ) {
   nObs <- nrow( data )
   coef <- aidsCoef( allCoef )
   nGoods <- length( coef$alpha )
   hom <- all.equal( rowSums( coef$gamma ), rep( 0, nGoods ) ) == TRUE
   sym <- all.equal( coef$gamma, t( coef$gamma ) ) == TRUE
   lnp <- aidsPx( "TL", pNames, coef = coef, data = data, alpha0 = alpha0 )
   result <- matrix( 0, nrow = nObs * ( nGoods - 1 ),
      ncol = ( nGoods + 2 ) * ( nGoods - 1 ) )
   for( eq in 1:( nGoods - 1 ) ) {
      myRows <- ( ( eq - 1 ) * nObs + 1 ):( eq * nObs )
      # derivatives of alphas
      for( i in 1:( nGoods - 1 ) ) {
         myCol <- ( i - 1 ) * ( nGoods + 2 ) + 1
         result[ myRows, myCol ] <- ( i == eq ) -
            coef$beta[ eq ] *
            ( log( data[[ pNames[ i ] ]] ) -
            log( data[[ pNames[ nGoods ] ]] ) )
      }
      # derivatives of betas
      myCol <- ( eq - 1 ) * ( nGoods + 2 ) + 2
      result[ myRows, myCol ] <- log( data[[ xtName ]] ) - lnp
      # derivatives of gammas
      for( i in 1:( nGoods - 1 ) ) {
         for( j in 1:( nGoods - hom ) ) {
            myCol <- ( i - 1 ) * ( nGoods + 2 ) + 2 + j
            result[ myRows, myCol ] <-
               ( i == eq ) * ( log( data[[ pNames[ j ] ]] ) -
                  hom * log( data[[ pNames[ nGoods ] ]] ) ) -
               0.5 * coef$beta[ eq ] *
               ( log( data[[ pNames[ i ] ]] ) -
                  log( data[[ pNames[ nGoods ] ]] ) ) *
               ( log( data[[ pNames[ j ] ]] ) -
                  hom * log( data[[ pNames[ nGoods ] ]] ) )
         }
      }
   }
   delCols <- NULL
   for( i in 1:( nGoods - 1 ) ) {
      if( hom ) {
         delCols <- c( delCols, i * ( nGoods + 2 ) )
      }
      if( sym && i >= 2 ) {
         for( j in 1:( i - 1 ) ) {
            delCol <- ( i - 1 ) * ( nGoods + 2 ) + 2 + j
            stayCol <- ( j - 1 ) * ( nGoods + 2 ) + 2 + i
            result[ , stayCol ] <- result[ , stayCol ] + result[ , delCol ]
            delCols <- c( delCols, delCol )
         }
      }
   }
   if( !is.null( delCols ) ) {
      result <- result[ , -delCols ]
   }
   return( result )
}
