aidsRestr <- function( nGoods, hom = TRUE, sym = TRUE, LA = TRUE ) {
   if( sym && !hom ) {
      hom <- TRUE  # symmetry implies homogeneity
      warning( "symmetry implies homogeneity: imposing additionally homogeniety" )
   }
   restr <- NULL
   if( LA ) {
      if( hom ) {
         restr <- matrix( 0, nGoods - 1, ( nGoods - 1 ) * ( nGoods + 2 ) )
         for( i in 1:( nGoods - 1 ) ) {
            for( j in 1:nGoods ) {
               restr[ i, ( i - 1 ) * ( nGoods + 2 ) + 2 + j ] <- 1
            }
         }
      }
      if( sym ) {
         restr <- rbind( restr, matrix( 0, ( nGoods - 1 ) * ( nGoods - 2 ) / 2,
            ( nGoods - 1 ) * ( nGoods + 2 ) ) )
         k <- 0
         for( i in 1:( nGoods - 2 ) ) {
            for( j in ( i + 1 ):( nGoods - 1 ) ) {
               k <- k + 1
               restr[ nGoods - 1 + k, ( i - 1 ) * ( nGoods + 2 ) +
                  2 + j ] <-  1
               restr[ nGoods - 1 + k, ( j - 1 ) * ( nGoods + 2 ) +
                  2 + i ] <- -1
            }
         }
      }
   }
   return( restr )
}
