snqProfitEst <- function( pNames, qNames, fNames = NULL,
   ivNames = NULL, data,  form = 0, base = 1,
   weights = snqProfitWeights( pNames, qNames, data, "DW92", base = base ),
   method = ifelse( is.null( ivNames ), "SUR", "3SLS" ), ... ) {

   checkNames( c( pNames, qNames, fNames, ivNames ), names( data ) )

   if( length( qNames ) != length( pNames ) ) {
      stop( "arguments 'qNames' and 'pNames' must have the same length." )
   }
   if( length( pNames ) < 2 ) {
      stop( "you must specify at least 2 netputs." )
   }
   if( length( pNames ) != length( weights ) ) {
      stop( "arguments 'pNames' and 'weights' must have the same length." )
   }
   if( min( weights ) < 0 ) {
      warning( paste( "At least one weight of the prices for normalization",
         "(argument 'weights') is negative. Thus, in this case positive",
         "semidefiniteness of the 'beta' matrix does not ensure",
         "a convex profit function." ) )
   }

   nNetput <- length( qNames )  # number of netputs
   nFix    <- length( fNames )  # number of fixed inputs
   nIV     <- length( ivNames )  # number of fixed inputs
   nObs    <- nrow( data )      # number of observations

   if( form == 0 ) {
      nCoef   <- nNetput + nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
         ( nFix + 1 ) * nFix/2  #number of coefficients
   } else if( form == 1 ) {
      nCoef   <- nNetput + nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
         nNetput * ( nFix + 1 ) * nFix/2  #number of coefficients
   } else {
      stop( "argument 'form' must be either 0 or 1." )
   }

   result  <- list()

   ## scaling data
   if( !is.null( base ) ) {
      estData <- data.frame( nr = 1:nObs )
      for( i in 1:nNetput ) {
         estData[[ pNames[ i ] ]] <- data[[ pNames[ i ] ]] /
            mean( data[[ pNames[ i ] ]][ base ] )
         estData[[ qNames[ i ] ]] <- data[[ qNames[ i ] ]] *
            mean( data[[ pNames[ i ] ]][ base ] )
      }
      if( !is.null( fNames ) ) {
         for( i in 1:nFix ) {
            estData[[ fNames[ i ] ]] <- data[[ fNames[ i ] ]]
         }
      }
   } else {
      estData <- data
   }

   ## price index for normalization
   modelData <- data.frame( nr = 1:nObs, normPrice = 0 )
   for( i in 1:nNetput ) {
      modelData$normPrice <- modelData$normPrice +
         estData[[ pNames[ i ] ]] * weights[ i ]
   }

   ## real/normalized netput prices and netput quantities
   result$pMeans <- array( NA, nNetput )
   result$qMeans <- array( NA, nNetput )
   for( i in 1:nNetput ) {
      modelData[[ paste( "pr", as.character( i ), sep = "" ) ]] <-
         estData[[ pNames[ i ] ]] / modelData$normPrice
      result$pMeans[ i ] <- mean( estData[[ pNames[ i ] ]] )
      modelData[[ paste( "q", as.character( i ), sep = "" ) ]] <-
         estData[[ qNames[ i ] ]]
      result$qMeans[ i ] <- mean( estData[[ qNames[ i ] ]] )
   }
   names( result$pMeans ) <- pNames
   names( result$qMeans ) <- qNames

   ## quadratic netput prices
   for( i in 1:nNetput ) {
      for( j in 1:nNetput ) {
         for( k in 1:nNetput ) {
            modelData[[ paste( "pq", as.character( i ), ".", as.character( j ), ".",
               as.character( k ), sep = "" ) ]] <-
               -0.5 * weights[ i ] * estData[[ pNames[ j ] ]] *
               estData[[ pNames[ k ] ]] / modelData$normPrice^2
         }
      }
   }

   ## quasi-fix inputs
   if( nFix > 0 ) {
      for( i in 1:nFix ) {
         result$fMeans[ i ] <- mean( data[[ fNames[ i ] ]] )
         modelData[[ paste( "f", as.character( i ), sep = "" ) ]] <-
            data[[ fNames[ i ] ]] / result$fMeans[ i ]
      }
      names( result$fMeans ) <- fNames
      ## quadratic quasi-fix inputs
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            for( k in 1:nFix ) {
               modelData[[ paste( "fq", as.character( i ), ".", as.character( j ), ".",
               as.character( k ), sep = "" ) ]] <-
                  0.5 * ifelse( form == 0, weights[ i ], 1 ) *
                  ( data[[ fNames[ j ] ]] / result$fMeans[ j ] ) *
                  ( data[[ fNames[ k ] ]] / result$fMeans[ k ] )
            }
         }
      }
   }

   ## instrumental variables
   if( nIV == 0 ) {
      inst <- NULL
   } else {
      inst <- as.formula( paste( "~", paste( ivNames, collapse = "+" ) ) )
      for( i in 1:nIV ) {
         modelData[[ paste( "iv", as.character( i ), sep = "" ) ]] <-
            data[[ ivNames[ i ] ]]
      }
   }
   system <- snqProfitSystem( nNetput, nFix )    # equation system
   restrict <- snqProfitRestrict( nNetput, nFix, form )    # restrictions
   result$est <- systemfit( method = method, eqns = system, data = modelData,
      TX = restrict, inst = inst, ... )
   result$coef <- snqProfitCoef( coef = result$est$bt, nNetput = nNetput,
      nFix = nFix, form = form, coefCov = result$est$btcov,
      df = nNetput * nObs - nCoef, 
      qNames = qNames, pNames = pNames, fNames = fNames )
      # estimated coefficients
   result$coef$liCoef <- result$est$bt
   result$coef$liCoefCov <- result$est$btcov

   if( nFix > 0 ) {
      # delta[ i, j ]
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            result$coef$delta[ i, j ] <- result$coef$delta[ i, j ] /
               result$fMeans[ j ]
            # all coefficients
            k <- nNetput + nNetput^2 + ( i - 1 ) * nFix + j
            result$coef$allCoef[ k ] <- result$coef$allCoef[ k ] /
               result$fMeans[ j ]
            result$coef$allCoefCov[ k, ] <- result$coef$allCoefCov[ k, ] /
               result$fMeans[ j ]
            result$coef$allCoefCov[ , k ] <- result$coef$allCoefCov[ , k ] /
               result$fMeans[ j ]
            result$coef$stats[ k, c( 1, 2 ) ] <- result$coef$stats[ k, c( 1, 2 ) ] /
               result$fMeans[ j ]
            # linear independent coefficients
            k <- nNetput + ( nNetput * ( nNetput - 1 ) ) / 2 +
               ( i - 1 ) * nFix + j
            result$coef$liCoef[ k ] <- result$coef$liCoef[ k ] /
               result$fMeans[ j ]
            result$coef$liCoefCov[ k, ] <- result$coef$liCoefCov[ k, ] /
               result$fMeans[ j ]
            result$coef$liCoefCov[ , k ] <- result$coef$liCoefCov[ , k ] /
               result$fMeans[ j ]
         }
      }
      if( form == 0 ) {
         # gamma[ i, j ]
         for( i in 1:nFix ) {
            for( j in 1:nFix ) {
               # delta[ i, j ]
               result$coef$gamma[ i, j ] <- result$coef$gamma[ i, j ] /
                  ( result$fMeans[ i ] * result$fMeans[ j ] )
               # all coefficients
               k <- nNetput + nNetput^2 + nNetput * nFix + ( i - 1 ) * nFix + j
               result$coef$allCoef[ k ] <- result$coef$allCoef[ k ] /
                  ( result$fMeans[ i ] * result$fMeans[ j ] )
               result$coef$allCoefCov[ k, ] <- result$coef$allCoefCov[ k, ] /
                  ( result$fMeans[ i ] * result$fMeans[ j ] )
               result$coef$allCoefCov[ , k ] <- result$coef$allCoefCov[ , k ] /
                  ( result$fMeans[ i ] * result$fMeans[ j ] )
               result$coef$stats[ k, c( 1, 2 ) ] <- result$coef$stats[ k, c( 1, 2 ) ] /
                  ( result$fMeans[ i ] * result$fMeans[ j ] )
               # linear independent coefficients
               if( j >= i ) {
                  k <- nNetput + nNetput * (nNetput - 1 ) / 2 +
                     nNetput * nFix + veclipos( i, j, nFix )
                  result$coef$liCoef[ k ] <- result$coef$liCoef[ k ] /
                     ( result$fMeans[ i ] * result$fMeans[ j ] )
                  result$coef$liCoefCov[ k, ] <- result$coef$liCoefCov[ k, ] /
                     ( result$fMeans[ i ] * result$fMeans[ j ] )
                  result$coef$liCoefCov[ , k ] <- result$coef$liCoefCov[ , k ] /
                     ( result$fMeans[ i ] * result$fMeans[ j ] )
               }
            }
         }
      } else {
         # gamma[ n, i, j ]
         for( n in 1:nNetput ) {
            for( i in 1:nFix ) {
               for( j in 1:nFix ) {
                  # delta[ i, j ]
                  result$coef$gamma[ n, i, j ] <- result$coef$gamma[ n, i, j ] /
                     ( result$fMeans[ i ] * result$fMeans[ j ] )
                  # all coefficients
                  k <- nNetput + nNetput^2 + nNetput * nFix +
                     ( n - 1 ) * nFix^2 + ( i - 1 ) * nFix + j
                  result$coef$allCoef[ k ] <- result$coef$allCoef[ k ] /
                     ( result$fMeans[ i ] * result$fMeans[ j ] )
                  result$coef$allCoefCov[ k, ] <- result$coef$allCoefCov[ k, ] /
                     ( result$fMeans[ i ] * result$fMeans[ j ] )
                  result$coef$allCoefCov[ , k ] <- result$coef$allCoefCov[ , k ] /
                     ( result$fMeans[ i ] * result$fMeans[ j ] )
                  result$coef$stats[ k, c( 1, 2 ) ] <- result$coef$stats[ k, c( 1, 2 ) ] /
                     ( result$fMeans[ i ] * result$fMeans[ j ] )
                  # linear independent coefficients
                  if( j >= i ) {
                     k <- nNetput + ( nNetput * (nNetput - 1 ) ) / 2 +
                        nNetput * nFix + ( n - 1 ) * nFix * ( nFix + 1 ) / 2 +
                        veclipos( i, j, nFix )
                     result$coef$liCoef[ k ] <- result$coef$liCoef[ k ] /
                        ( result$fMeans[ i ] * result$fMeans[ j ] )
                     result$coef$liCoefCov[ k, ] <- result$coef$liCoefCov[ k, ] /
                        ( result$fMeans[ i ] * result$fMeans[ j ] )
                     result$coef$liCoefCov[ , k ] <- result$coef$liCoefCov[ , k ] /
                        ( result$fMeans[ i ] * result$fMeans[ j ] )
                  }
               }
            }
         }
      }
   }

   result$fitted <- data.frame( profit0 = rep( 0, nObs ) )
   for( i in 1:nNetput ) {
      result$fitted[[ qNames[ i ] ]] <- result$est$eq[[ i ]]$fitted
      result$fitted[[ "profit0" ]] <- result$fitted[[ "profit0" ]] +
         result$fitted[[ qNames[ i ] ]] * estData[[ pNames[ i ] ]]
   }
   result$fitted[[ "profit" ]] <- result$fitted[[ "profit0" ]]
   result$fitted[[ "profit0" ]] <- NULL

   result$r2 <- array( NA, c( nNetput ) )
   for( i in 1:nNetput ) {
      result$r2[ i ] <- result$est$eq[[ i ]]$r2
   }
   names( result$r2 ) <- qNames

   result$hessian <- snqProfitHessian( result$coef$beta, result$pMeans, weights )
      # Hessian matrix
   result$ela <- snqProfitEla( result$coef$beta, result$pMeans,
      result$qMeans, weights )   # estimated elasticities
   result$estData  <- estData
   result$weights  <- weights
   result$normPrice <- modelData$normPrice
   result$convexity  <- semidefiniteness( result$hessian[
      1:( nNetput - 1 ), 1:( nNetput - 1 ) ] )$positive
   result$form <- form
   class( result )  <- "snqProfitEst"
   return( result )
}
