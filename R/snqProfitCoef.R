snqProfitCoef <- function( coef, nNetput, nFix, form = 0, coefCov = NULL, 
   df = 1, pNames = NULL, qNames = NULL, fNames = NULL ) {

   nCoef    <- nNetput + nNetput * ( nNetput - 1 )/2 + nNetput * nFix
      # number of coefficients
   if( form == 0 ) {
      nCoef <- nCoef + ( nFix + 1 ) * nFix / 2
   } else if ( form == 1 ) {
      nCoef <- nCoef + nNetput * ( nFix + 1 ) * nFix / 2
   } else {
      stop( "argument 'form' must be either 0 or 1" )
   }
   if( length( coef ) != nCoef ) {
      stop( paste( "A SNQ profit function of form =", form, "with", nNetput,
         "Netputs and", nFix, "fix inputs must have exactly", nCoef,
         "linear independent coefficients (argument 'coef')." ) )
   }

   result <- list()
   result$alpha <- coef[ 1:nNetput ]
   result$beta <- array( 0, c( nNetput, nNetput ) )
   for( i in 1:( nNetput - 1 ) ) {
      for( j in 1:( nNetput - 1 ) ) {
         result$beta[ i, j ] <- coef[ nNetput + veclipos( i, j, nNetput - 1 ) ]
         result$beta[ i, nNetput ] <- result$beta[ i, nNetput ] - result$beta[ i, j ]
      }
      result$beta[ nNetput, i ] <- result$beta[ i, nNetput ]
      result$beta[ nNetput, nNetput ] <- result$beta[ nNetput, nNetput ] -
         result$beta[ i, nNetput ]
   }
   if( nFix > 0 ) {
      result$delta <- array( 0, c( nNetput, nFix ) )
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            result$delta[ i, j ] <- coef[ nNetput + nNetput * ( nNetput - 1 ) / 2 +
               ( i - 1 ) * nFix + j ]
         }
      }
      if( form == 0 ) {
         result$gamma <- array( 0, c( nFix, nFix ) )
         for( i in 1:nFix ) {
            for( j in 1:nFix ) {
               result$gamma[i,j] <- coef[ nNetput + nNetput * ( nNetput - 1 ) / 2 +
                  nNetput * nFix + veclipos( i, j, nFix ) ]
            }
         }
      } else if( form == 1 ) {
         result$gamma <- array( 0, c( nNetput, nFix, nFix ) )
         for( i in 1:nNetput ) {
            for( j in 1:nFix ) {
               for( k in 1:nFix ) {
                  result$gamma[ i, j, k ] <- coef[ nNetput +
                     nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
                     ( i - 1 ) * ( nFix + 1 ) * nFix / 2 +
                     veclipos( j, k, nFix ) ]
               }
            }
         }
      } else {
         stop( "argument 'form' must be either 0 or 1." )
      }
      result$allCoef <- c( result$alpha, array( result$beta ),
         array( t( result$delta ) ) )
      if( form == 0 ) {
         result$allCoef <- c( result$allCoef, array( result$gamma ) )
      } else {
         for( i in 1:nNetput ) {
            result$allCoef <- c( result$allCoef, array( result$gamma[ i, , ] ) )
         }
      }
   } else {
      result$allCoef <- c( result$alpha, array( result$beta ) )
   }
   names( result$allCoef ) <- snqProfitCoefNames( nNetput, nFix,
         form = form, all = TRUE )
   ## completing the coefficient covariance matrix
   if( !is.null( coefCov ) ) {
      result$allCoefCov <- coefCov  # covariance matrix of all coefficients
      for( i in 1:nNetput ) {
         if( i >= 2 ) {
            for( j in 1:( i - 1 ) ) {
               # beta(j,i)
               k  <- i * nNetput + j
               k2 <- j * nNetput + i
               result$allCoefCov <- insertRow( result$allCoefCov, k )
               result$allCoefCov[ k, ] <- result$allCoefCov[ k2, ]
               result$allCoefCov <- insertCol( result$allCoefCov, k )
               result$allCoefCov[ , k ] <- result$allCoefCov[ , k2 ]
            }
         }
         # beta(i,nNetput)
         k  <- ( 1 + i ) * nNetput
         k2 <- i * nNetput + 1
         result$allCoefCov <- insertRow( result$allCoefCov, k )
         if( nNetput > 2 ) {
            result$allCoefCov[ k, ] <- -colSums( result$allCoefCov[
               k2:( k - 1 ), ] )
         } else {
            result$allCoefCov[ k, ] <- -result$allCoefCov[ k2, ]
         }
         result$allCoefCov <- insertCol( result$allCoefCov, k )
         if( nNetput > 2 ) {
            result$allCoefCov[ , k ] <- -rowSums( result$allCoefCov[ ,
               k2:( k - 1 ) ] )
         } else {
            result$allCoefCov[ , k ] <- -result$allCoefCov[ , k2 ]
         }
      }
      if( nFix > 0 ) {
         if( form == 0 ) {
            for( i in 1:nFix ) {
               if( i >= 2 ) {
                  for( j in 1:(i-1) ) {
                     # gamma(j,i)
                     k  <- nNetput + nNetput^2 + nNetput * nFix + (i-1) * nFix + j
                     k2 <- nNetput + nNetput^2 + nNetput * nFix + (j-1) * nFix + i
                     result$allCoefCov <- insertRow( result$allCoefCov, k )
                     result$allCoefCov[ k, ] <- result$allCoefCov[ k2, ]
                     result$allCoefCov <- insertCol( result$allCoefCov, k )
                     result$allCoefCov[ , k] <- result$allCoefCov[ ,k2 ]
                  }
               }
            }
         } else if( form == 1 ) {
            for( n in 1:nNetput ) {
               for( i in 1:nFix ) {
                  if( i >= 2 ) {
                     for( j in 1:(i-1) ) {
                        # gamma(j,i)
                        k  <- nNetput + nNetput^2 + nNetput * nFix +
                           ( n - 1 ) * nFix^2 + ( i - 1 ) * nFix + j
                        k2 <- nNetput + nNetput^2 + nNetput * nFix +
                           ( n - 1 ) * nFix^2 + ( j - 1 ) * nFix + i
                        result$allCoefCov <- insertRow( result$allCoefCov, k )
                        result$allCoefCov[ k, ] <- result$allCoefCov[ k2, ]
                        result$allCoefCov <- insertCol( result$allCoefCov, k )
                        result$allCoefCov[ , k] <- result$allCoefCov[ ,k2 ]
                     }
                  }
               }
            }
         } else {
            stop( "argument 'form' must be either 0 or 1." )
         }
      }
      rownames( result$allCoefCov ) <- snqProfitCoefNames( nNetput, nFix,
         form = form, all = TRUE )
      colnames( result$allCoefCov ) <- snqProfitCoefNames( nNetput, nFix,
         form = form, all = TRUE )
      result$stats <- array( NA, c( length( result$allCoef ), 4 ) )
      result$stats[ , 1 ] <- result$allCoef
      result$stats[ , 2 ] <- sqrt( diag( result$allCoefCov ) )
      result$stats[ , 3 ] <- result$stats[ , 1 ] / result$stats[ , 2 ]
      result$stats[ , 4 ] <- 2 * ( 1 - pt( abs( result$stats[ , 3 ] ), df ) )
      rownames( result$stats ) <- snqProfitCoefNames( nNetput, nFix,
         form = form, all = TRUE )
      colnames( result$stats ) <- c( "value", "std.err", "t-value", "prob" )
   }
   if( !is.null( qNames ) ) {
      if( length( qNames ) != nNetput ) {
         stop( paste( "argument 'qNames' must have as many elements as",
            "there are netputs" ) )
      }
      names( result$alpha ) <- qNames
      rownames( result$beta ) <- qNames
      if( nFix > 0 ) {
         rownames( result$delta ) <- qNames
      }
   }
   if( !is.null( pNames ) ) {
      if( length( pNames ) != nNetput ) {
         stop( paste( "argument 'pNames' must have as many elements as",
            "there are netputs" ) )
      }
      colnames( result$beta ) <- pNames
   }
   if( !is.null( fNames ) ) {
      if( length( fNames ) != nFix ) {
         stop( paste( "argument 'fNames' must have as many elements as",
            "there are fixed inputs" ) )
      }
      colnames( result$delta ) <- fNames
      if( form == 0 ) {
         rownames( result$gamma ) <- fNames
         colnames( result$gamma ) <- fNames
      }
   }
   return( result )
}
