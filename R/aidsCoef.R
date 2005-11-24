aidsCoef <- function( b, cov = NULL, df = 1, LA = TRUE,
   pNames = NULL, wNames = NULL ) {
   nGoods <- -0.5 + ( 2.25 + nrow( array( b ) ) )^0.5
   if( LA ) {
      M <- matrix( 0, nGoods * ( nGoods + 2 ), ( nGoods - 1 ) * ( nGoods + 2 ) )
      for(i in 1:( nGoods - 1 ) ) {
         M[ i,( i - 1 ) * ( nGoods + 2 ) + 1 ]   <-  1   # alpha[i]
         M[ nGoods, ( i - 1 ) * ( nGoods + 2 ) + 1 ]   <- -1   # alpha[nGoods]
         M[ nGoods + i, ( i - 1 ) * ( nGoods + 2 ) + 2] <-  1   # beta[i]
         M[ nGoods + nGoods, ( i - 1 ) * ( nGoods + 2 ) + 2 ] <- -1
            # beta[ nGoods ]
         for( j in 1:nGoods ) {
            M[ 2 * nGoods + ( i - 1 ) * nGoods + j,
               (i - 1 ) * ( nGoods + 2 ) + 2 + j ] <-  1   # gamma[i,j]
            M[ 2 * nGoods + ( nGoods - 1 ) * nGoods + j,
               ( i - 1 ) * ( nGoods + 2 ) + 2 + j ] <- -1   # gamma[nGoods,j]
         }
      }
      all     <- M %*% b
      all[nGoods]  <- all[nGoods]+1
      names( all ) <- c(
            paste( "alpha", c( 1:nGoods ) ),
            paste( "beta", c( 1:nGoods ) ),
            paste( "gamma", rep( 1:nGoods, each = nGoods ),
               rep( 1:nGoods, nGoods ) ) )
      alpha   <- all[1:nGoods]
      beta    <- all[(nGoods+1):(2*nGoods)]
      gamma   <- t(array(all[(2*nGoods+1):(nGoods*(nGoods+2))],c(nGoods,nGoods)))
      allcov <- NULL
      stat    <- NULL
      if(!is.null(cov)) {
         allcov   <- M %*% cov %*% t(M)
         stat     <- coefTable( all, sqrt( diag( allcov ) ), df )
      }
   }
   if( !is.null( wNames ) ) {
      names( alpha ) <- wNames
      names( beta ) <- wNames
      rownames( gamma ) <- wNames
   }
   if( !is.null( pNames ) ) {
      colnames( gamma ) <- pNames
   }
   coef <- list()
   coef$alpha <- alpha
   coef$beta  <- beta
   coef$gamma <- gamma
   coef$all   <- all
   coef$allcov <- allcov
   coef$stat  <- stat
   return( coef )
}
