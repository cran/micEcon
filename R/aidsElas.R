aidsElas <- function( coef, shares, prices = NULL, method = "AIDS",
   quantNames = NULL, priceNames = NULL, coefVcov = NULL, df = NULL ) {

   nGoods <- length( coef$alpha )

   if( length( coef$alpha ) != length( coef$beta ) ) {
      stop( "arguments 'alpha' and 'beta' must have the same length" )
   } else if( nrow( coef$gamma ) != ncol( coef$gamma ) ) {
      stop( "argument 'gamma' must be a square matrix" )
   } else if( length( coef$alpha ) != nrow( coef$gamma ) ) {
      stop( "number of rows of argument 'gamma' must be equal",
         " to the length of argument 'alpha'" )
   } else if(  length( coef$alpha ) != length( shares ) ) {
      stop( "arguments 'alpha' and 'shares' must have the same length" )
   } else if(  length( coef$alpha ) != length( prices ) && !is.null( prices ) ) {
      stop( "arguments 'alpha' and 'prices' must have the same length" )
   }
   if( is.null( quantNames ) ) {
      quantNames <- .aidsQuantNames( shares, coef, nGoods )
   } else {
      if( length( quantNames ) != nGoods ) {
         stop( "argument 'quantNames' must have ", nGoods, " elements" )
      }
   }
   if( is.null( priceNames ) ) {
      priceNames <- .aidsPriceNames( prices, coef, nGoods )
   } else {
      if( length( priceNames ) != nGoods ) {
         stop( "argument 'priceNames' must have ", nGoods, " elements" )
      }
   }


   if( method %in% c( "AIDS", "GA", "B1", "B2" ) ) {
      if( is.null( prices ) ) {
         stop( "methods 'AIDS', 'GA', 'B1', and 'B2'",
            " require argument 'prices'" )
      }
   } else if( method %in% c( "Go", "Ch", "EU" ) ) {
      if( !is.null( prices ) ) {
         warning( "methods 'Go', 'Ch', and 'EU'",
            " do not require argument 'prices'" )
      }
   }

   ela <- list()
   ela$method <- method

   ones <- rep( 1, nGoods )

   if( method == "AIDS" ) {
      ela$exp <- ones + coef$beta/shares
      ela$hicks <- -diag( 1, nGoods, nGoods ) +
         ones %*% t( shares ) +
         coef$gamma / ( shares %*% t( ones ) ) -
         coef$beta %*% t( ones ) *
         ( ones %*% t( coef$alpha ) -
         ones %*% t( shares )+
         ones %*% t( coef$gamma %*% log( prices ))) /
         ( shares %*% t( ones ) )
      ela$marshall <- ela$hicks - ( ela$exp %*% t( ones ) ) *
         ( ones %*% t( shares ))
   } else if( method %in% c( "Ch", "Go" ) ) {
      ela$exp <- ones + coef$beta / shares
      ela$hicks <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( shares %*% t( ones ) ) +
         ones %*% t( shares )
      ela$marshall <- ela$hicks - ( ela$exp %*% t( ones ) ) *
         ( ones %*% t(shares))
   } else if( method == "EU" ) {
      ela$exp <- ones + coef$beta / shares
      ela$marshall <- -diag( 1, nGoods, nGoods ) + coef$gamma /
         ( shares %*% t( ones ) )
      ela$hicks <- ela$marshall + ( ela$exp %*% t( ones ) ) *
         ( ones %*% t(shares))
   } else if( method %in% c( "GA", "B1" ) ) {
      denom <- 1  # part of denominator for exp. + Marsh. elasticities
      for( k in 1:nGoods ) {
         denom <- denom + coef$beta[ k ] * log( prices[ k ] )
      }
      ela$exp <- ones + coef$beta / ( shares * ( 1 + denom ) )
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            numer <- shares[ j ] # part of numerator for Marsh. elasticities
            for( k in 1:nGoods ) {
               numer <- numer + coef$gamma[ k, j ] * log( prices[ k ] )
            }
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               ( coef$beta[ i ] / shares[ i ] ) * numer / denom
         }
      }
      ela$hicks <- ela$marshall + ( ela$exp %*% t( ones ) ) *
         ( ones %*% t(shares))
   } else if( method %in% c( "B2" ) ) {
      paren <- 1  # term in parenthesis for expenditure elasticities
      for( k in 1:nGoods ) {
         paren <- paren - coef$beta[ k ] * log( prices[ k ] )
      }
      ela$exp <- ones + ( coef$beta / shares ) * paren
      ela$marshall <- matrix( NA, nGoods, nGoods )
      for( i in 1:nGoods ) {
         for( j in 1:nGoods ) {
            parenSmall <- coef$alpha[ j ] # term in small par. for Marsh. elast.
            for( l in 1:nGoods ) {
               parenSmall <- parenSmall +
                  coef$gamma[ l, j ] * log( prices[ l ] )
            }
            parenBig <- shares[ j ] # term in big parenthesis for Marsh. elast.
            for( k in 1:nGoods ) {
               parenBig <- parenBig +
                  coef$gamma[ k, j ] * log( prices[ k ] ) -
                  coef$beta[ k ] * log( prices[ k ] ) * parenSmall
            }
            ela$marshall[ i, j ] <- -( i == j ) + coef$gamma[ i, j ] / shares[ i ] -
               ( coef$beta[ i ] / shares[ i ] ) * parenBig
         }
      }
      ela$hicks <- ela$marshall + ( ela$exp %*% t( ones ) ) *
         ( ones %*% t(shares))
   } else {
      stop( "argument 'method' must be either 'AIDS', 'GA', 'B1', 'B2',",
         " 'Go', 'Ch', or 'EU'" )
   }
   names( ela$exp )         <- quantNames
   rownames( ela$hicks )    <- quantNames
   colnames( ela$hicks )    <- priceNames
   rownames( ela$marshall ) <- quantNames
   colnames( ela$marshall ) <- priceNames
   if( !is.null( coefVcov ) && method %in% c( "AIDS" ) ) {
      jacobian <- aidsElasJacobian( coef = coef, share = shares, price = prices,
         method = method, quantNames = quantNames, priceNames = priceNames )
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
      ela$df <- df
   }
   class( ela ) <- "aidsElas"
   return( ela )
}
