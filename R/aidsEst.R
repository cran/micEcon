aidsEst <- function( pNames, wNames, xtName,
      data = NULL, ivNames = NULL, qNames = wNames,
      method = "LA:L", hom = TRUE, sym = TRUE,
      elaFormula = "Ch", pxBase = 1,
      estMethod = ifelse( is.null( ivNames ), "SUR", "3SLS" ),
      maxiterIL = 50, tolIL = 1e-5, alpha0 = 0, TX = FALSE, ... ) {

   if( length( pNames ) != length( wNames ) ) {
      stop( "arguments 'pNames' and 'wNames' must have the same length" )
   }
   nGoods <- length( pNames )
   extractPx <- function( method ) {
      px <- substr( method, 4, nchar( method ) )
      if( !( px %in% c( "S", "SL", "P", "L", "T" ) ) ) {
         stop( "no valid price index specified!" )
      }
      return( px )
   }

   if( substr( method, 1, 2 ) == "LA" ) {
      if( nchar( method ) < 4 ) {
         warning( "No price index specified: using Laspeyres price index" )
         px <- "L"
      } else {
         px <- extractPx( method )
      }
   } else if ( substr( method, 1, 2 ) %in% c( "MK", "IL" ) ) {
      if( nchar( method ) < 4 ) {
         warning( "No initial price index specified:",
            " using Laspeyres price index" )
         px <- "L"
      } else {
         px <- extractPx( method )
      }
   } else {
      stop( "at the moment only the methods",
         " 'Linear Approximation' (LA) and",
         " 'Iterated Linear Least Squares' (IL)",
         " are supported" )
   }
   if( sym && !hom ) {
      hom <- TRUE  # symmetry implies homogeneity
      warning( "symmetry implies homogeneity: imposing additionally homogeniety" )
   }
   allVarNames <- c( pNames, wNames, xtName, ivNames )
   if( sum( is.na( data[ , allVarNames ] ) ) > 0 ) {
      warning( "there are some NAs in the data,",
         " all observations (rows) with NAs are excluded from the analysis" )
      data <- data[ !is.na( rowSums( data[ , allVarNames ] ) ), ]
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
   sysData <- data.frame( xt = data[[ xtName ]],
      lxtr = ( log( data[[ xtName ]] ) - lnp ) )
   for( i in 1:nGoods ) {
      sysData[[ paste( "w", i, sep = "" ) ]] <- data[[ wNames[ i ] ]]
      sysData[[ paste( "lp", i, sep = "" ) ]] <- log( data[[ pNames[ i ] ]] )
   }
   if( is.null( ivNames )) {
      ivFormula <- NULL
   } else {
      estMethod <- "3SLS"
      ivFormula <- "~"
      for( i in 1:length( ivNames ) ) {
         sysData[[ paste( "i", i, sep = "" ) ]] <- data[[ ivNames[ i ] ]]
         ivFormula <- paste( ivFormula, " + i", as.character( i ), sep = "" )
      }
      ivFormula <- as.formula( ivFormula )
   }
   restr <- aidsRestr( nGoods, hom, sym, TX = TX )
      # restrictions for homogeneity and symmetry
   system <- aidsSystem( nGoods )    # LA-AIDS equation system
   # estimate system
   if( TX ) {
      est <- systemfit( estMethod, system, data = sysData, TX = restr,
         inst = ivFormula, ... )
   } else {
      est <- systemfit( estMethod, system, data = sysData, R.restr = restr,
         inst = ivFormula, ... )
   }
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
   } else if( substr( method, 1, 2 ) %in% c( "MK", "IL" ) ) {
      b       <- est$b      # coefficients
      bd      <- est$b      # difference of coefficients between
                            # this and previous step
      iter    <- est$iter   # iterations of each SUR estimation
      iterIL <- 1          # iterations of IL Loop
      while( ( ( t( bd ) %*% bd ) / ( t( b ) %*% b ) )^0.5 > tolIL &&
            iterIL < maxiterIL ) {
         iterIL <- iterIL + 1      # iterations of IL Loop
         bl     <- b              # coefficients of previous step
         sysData$lxtr <- log( data[[ xtName ]] ) -
            aidsPx( "TL", pNames, wNames, data = data,
            alpha0 = alpha0, coef = aidsCoef( est$b ) )
            # real total expenditure using Translog price index
         if( TX ) {
            est <- systemfit( estMethod, system, data = sysData, TX = restr,
               inst = ivFormula, ... )    # estimate system
         } else {
            est <- systemfit( estMethod, system, data = sysData, R.restr = restr,
               inst = ivFormula, ... )    # estimate system
         }
         iter <- c( iter, est$iter ) # iterations of each estimation
         b    <- est$b   # coefficients
         bd   <- b - bl  # difference between coefficients from this
                         # and previous step
      }
      # calculating log of "real" (deflated) total expenditure
      sysData$lxtr <- log( data[[ xtName ]] ) -
         aidsPx( "TL", pNames, data = data,
         alpha0 = alpha0, coef = aidsCoef( est$b ) )
      # calculating matrix G
      Gmat <- cbind( rep( 1, nObs ), sysData$lxtr )
      for( i in 1:( nGoods ) ) {
         Gmat <- cbind( Gmat, sysData[[ paste( "lp", i, sep = "" ) ]] )
      }
      # testing matrix G
      if( FALSE ) {
         for( i in 1:( nGoods - 1 ) ) {
            print( est$eq[[ i ]]$fitted - Gmat %*%
               est$b[ ( ( i - 1 ) * ( nGoods + 2 ) + 1 ):( i * (nGoods + 2 ) ) ] )
         }
      }
      # calculating matrix J
      jacobian <- aidsJacobian( est$b, pNames, xtName, data = data,
         alpha0 = alpha0 )
      if( hom ) {
         TXmat <- aidsRestr( nGoods, hom, sym, TX = TRUE )
      } else {
         TXmat <- diag( ( nGoods - 1 ) * ( nGoods + 2 ) )
      }
      # Jmat <- t( TXmat ) %*% ( diag( nGoods - 1 ) %x% t( Gmat ) ) %*% jacobian
      # JmatInv <- TXmat %*% solve( Jmat ) %*% t( TXmat )
      # bcov <- JmatInv  %*% ( est$rcov %x% ( t( Gmat ) %*% Gmat ) ) %*%
      #    t( JmatInv )
      Jmat <- crossprod( TXmat, ( diag( nGoods - 1 ) %x% t( Gmat ) ) ) %*% jacobian
      JmatInv <- TXmat %*% solve( Jmat, t( TXmat ) )
      bcov <- JmatInv  %*% ( est$rcov %x% crossprod( Gmat ) ) %*%
         t( JmatInv )
      result$coef <- aidsCoef( est$b, bcov, pNames = pNames,
         wNames = wNames, df = est$df )  # coefficients
      result$coef$alpha0 <- alpha0
      result$ela  <- aidsEla( result$coef, wMeans, pMeans,
         formula = "AIDS", pNames = pNames, qNames = qNames,
         coefVcov = result$coef$allcov, df = est$df )   # elasticities
      result$wFitted <- aidsCalc( pNames, xtName, data = data,
         coef = result$coef, alpha0 = alpha0, px = "TL" )$shares
         # estimated budget shares
      result$iterIL <- iterIL
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
   result$wMeans <- wMeans
   result$pMeans <- pMeans
   class( result ) <- "aidsEst"
   return( result )
}
