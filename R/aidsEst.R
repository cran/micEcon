aidsEst <- function( pNames, wNames, xtName,
      data = NULL, ivNames = NULL, qNames = wNames,
      method = "LA:L", hom = TRUE, sym = TRUE,
      elaFormula = "Ch", pxBase = 1,
      estMethod = ifelse( is.null( ivNames ), "SUR", "3SLS" ),
      maxiterMk = 50, tolMk = 1e-5, alpha0 = 0, ... ) {

   if( length( pNames ) != length( wNames ) ) {
      stop( "arguments 'pNames' and 'wNames' must have the same length." )
   }
   nGoods <- length( pNames )
   extractPx <- function( method ) {
      px <- substr( method, 4, nchar( method ) )
      if( !( px %in% c( "S", "SL", "P", "L", "T" ) ) ) {
         stop( "No valid price index specified!" )
      }
      return( px )
   }

   if( substr( method, 1, 2 ) != "LA" ) {
      if( nchar( method ) < 4 ) {
         warning( "No price index specified: using Laspeyres price index" )
         px <- "L"
      } else {
         px <- extractPx( method )
      }
   } else if ( substr( method, 1, 2 ) != "MK" ) {
      if( nchar( method ) < 4 ) {
         warning( paste( "No initial price index specified:",
            "using Laspeyres price index" ) )
         px <- "L"
      } else {
         px <- extractPx( method )
      }
   } else {
      stop( paste( "At the moment only the methods 'Linear Approximation'",
         "(LA) and 'Michalek & Keyzer' (MK) are supported" ) )
   }
   if( sym && !hom ) {
      hom <- TRUE  # symmetry implies homogeneity
      warning( "symmetry implies homogeneity: imposing additionally homogeniety" )
   }
   nObs   <- nrow( data )      # number of observations
   sample <- if( px == "SL") c( 2:nObs ) else c( 1:nObs )
   result <- list()
   wMeans <- numeric( nGoods )  # mean expenditure shares
   pMeans <- numeric( nGoods )  # mean prices
   for( i in seq( nGoods ) ) {
      wMeans[ i ] <- mean( data[[ wNames[ i ] ]][ sample ] )
      pMeans[ i ] <- mean( data[[ pNames[ i ] ]][ sample ] )
   }
   # log of price index
   lnp  <- aidsPx( px, pNames, wNames, data, base = pxBase )
   # prepare data.frame
   sysData <- data.frame( lxtr = ( log( data[[ xtName ]] ) - lnp ) )
   for( i in 1:nGoods ) {
      sysData <- cbind( sysData, data[[ wNames[ i ] ]] )
      names( sysData )[ length( sysData ) ] <-
         paste( "w", as.character( i ), sep = "" )

      sysData <- cbind( sysData, log( data[[ pNames[ i ] ]] ) )
      names( sysData )[ length( sysData ) ] <-
         paste( "lp", as.character( i ), sep = "" )
   }
   if( is.null( ivNames )) {
      ivFormula <- NULL
   } else {
      estMethod <- "3SLS"
      ivFormula <- "~"
      for( i in 1:length( ivNames ) ) {
         sysData <- cbind( sysData, data[[ ivNames[ i ] ]] )
         names( sysData )[ length( sysData ) ] <-
            paste( "i", as.character( i ), sep = "" )
         ivFormula <- paste( ivFormula, " + i", as.character( i ), sep = "" )
      }
      ivFormula <- as.formula( ivFormula )
   }
   restr <- aidsRestr( nGoods, hom, sym )
      # restrictions for homogeneity and symmetry
   system <- aidsSystem( nGoods )    # LA-AIDS equation system
   est <- systemfit( estMethod, system, data = sysData, R.restr = restr,
      inst = ivFormula, ... )   # estimate system
   if( substr( method, 1, 2 ) == "LA" ) {
      result$coef <- aidsCoef( est$b, est$bcov, pNames = pNames,
         wNames = wNames, df = est$df )   # coefficients
      if( !( elaFormula %in% c( "AIDS" ) ) ) {
         pMeans <- NULL
      }
      result$ela  <- aidsEla( result$coef, wMeans, pMeans,
         formula = elaFormula, pNames = pNames, qNames = qNames ) # elasticities
      result$wFitted <- aidsCalc( pNames, xtName, data = data,
         coef = result$coef, lnp = lnp )$shares   # estimated budget shares
      iter <- est$iter
   } else if( substr( method, 1, 2 ) == "MK" ) {
      b       <- est$b      # coefficients
      bd      <- est$b      # difference of coefficients between
                            # this and previous step
      iter    <- est$iter   # iterations of each SUR estimation
      iterMk <- 1          # iterations of M+K Loop
      while( ( ( t( bd ) %*% bd ) / ( t( b ) %*% b ) )^0.5 > tolMk &&
            iterMk < maxiterMk ) {
         iterMk <- iterMk + 1      # iterations of M+K Loop
         bl     <- b              # coefficients of previous step
         sysData$lxtr <- log( data[[ xtName ]] ) -
            aidsPx( "TL", pNames, wNames, data = data,
            alpha0 = alpha0, coef = aidsCoef( est$b ) )
            # real total expenditure using Translog price index
         est <- systemfit( estMethod, system, data = sysData, R.restr = restr,
            inst = ivFormula, ... )    # estimate system
         iter <- c( iter, est$iter ) # iterations of each estimation
         b    <- est$b   # coefficients
         bd   <- b - bl  # difference between coefficients from this
                         # and previous step
      }
      result$coef <- aidsCoef( est$b, est$bcov, pNames = pNames,
         wNames = wNames, df = est$df )  # coefficients
      result$coef$alpha0 <- alpha0
      result$ela  <- aidsEla( result$coef, wMeans, pMeans,
         formula = "AIDS", pNames = pNames, qNames = qNames )   # elasticities
      result$wFitted <- aidsCalc( pNames, xtName, data = data,
         coef = result$coef, alpha0 = alpha0, px = "TL" )$shares
         # estimated budget shares
      result$iterMk <- iterMk
   }
   names( result$wFitted ) <- paste( "wFitted", as.character( 1:nGoods ),
      sep = "" )
   result$wResid <- data.frame( matrix( NA, nrow = nObs, ncol = nGoods ) )
      # residuals of shares
   names( result$wResid ) <- paste( "wResid", as.character( 1:nGoods ), sep = "" )
   result$qObs <- data.frame( matrix( NA, nrow = nObs, ncol = nGoods ) )
      # observed quantities
   names( result$qObs ) <- paste( "qObs", as.character( 1:nGoods ), sep = "" )
   result$qFitted <- data.frame( matrix( NA, nrow = nObs, ncol = nGoods ) )
      # observed quantities
   names( result$qFitted ) <- paste( "qFitted", as.character( 1:nGoods ),
      sep = "" )
   result$qResid <- data.frame( matrix( NA, nrow = nObs, ncol = nGoods ) )
      # observed quantities
   names( result$qResid ) <- paste( "qResid", as.character( 1:nGoods ), sep = "" )
   for( i in 1:nGoods ) {
      result$wResid[ , i ] <- data[[ wNames[ i ] ]] - result$wFitted[ , i ]
      result$qObs[ , i ]   <- data[[ wNames[ i ] ]] * data[[ xtName ]] /
         data[[ pNames[ i ] ]]
      result$qFitted[ , i ] <- result$wFitted[ i ] * data[[ xtName ]] /
         data[[ pNames[ i ] ]]
      result$qResid[ , i ] <- result$qObs[ , i ] - result$qFitted[ , i ]
   }
   result$r2 <- array( 0, c( nGoods ) )
   for( i in 1:( nGoods - 1 ) ) {
      result$r2[ i ] <- est$eq[[ i ]]$r2
   }
   result$r2[ nGoods ] <- rSquared( data[[ wNames[ nGoods ] ]],
      result$wResid[ , nGoods ] )
   names( result$r2 ) <- wNames
   result$r2q <- array( 0, c( nGoods ) ) # R2 values for consumed quantities
   for( i in 1:nGoods ) {
      result$r2q[ i ] <- rSquared( result$qObs[ , i ], result$qResid[ , i ] )
   }
   result$iter <- iter
   result$est <- est
   result$method <- method
   result$px  <- px
   result$lnp <- lnp
   class( result ) <- "aidsEst"
   return( result )
}
