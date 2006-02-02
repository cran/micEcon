snqProfitEst <- function( pNames, qNames, fNames = NULL,
   ivNames = NULL, data,  form = 0, base = 1, scalingFactors = NULL,
   weights = snqProfitWeights( pNames, qNames, data, "DW92", base = base ),
   method = ifelse( is.null( ivNames ), "SUR", "3SLS" ), ... ) {

   checkNames( c( pNames, qNames, fNames, ivNames ), names( data ) )

   if( length( qNames ) != length( pNames ) ) {
      stop( "arguments 'qNames' and 'pNames' must have the same length" )
   }
   if( length( pNames ) < 2 ) {
      stop( "you must specify at least 2 netputs" )
   }
   if( length( pNames ) != length( weights ) ) {
      stop( "arguments 'pNames' and 'weights' must have the same length" )
   }
   if( min( weights ) < 0 ) {
      warning( "At least one weight of the prices for normalization",
         " (argument 'weights') is negative. Thus, in this case positive",
         " semidefiniteness of the 'beta' matrix does not ensure",
         " a convex profit function." )
   }
   if( !is.null( scalingFactors ) ) {
      if( length( scalingFactors ) != length( pNames ) ) {
         stop( "arguments 'pNames' and 'scalingFactors' must have the",
            " same length" )
      }
      if( base != 1 ) {
         warning( "argument 'base' is ignored because argument",
            " 'scalingFactors' is provided" )
      }
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
      stop( "argument 'form' must be either 0 or 1" )
   }

   result  <- list()

   ## scaling factors
   if( is.null( scalingFactors ) ) {
      scalingFactors <- rep( 1, nNetput )
      if( !is.null( base ) ) {
         for( i in 1:nNetput ) {
            scalingFactors[ i ] <- 1 / mean( data[[ pNames[ i ] ]][ base ] )
         }
      }
   }

   ## mean Values
   result$pMeans <- array( NA, nNetput )
   result$qMeans <- array( NA, nNetput )
   for( i in 1:nNetput ) {
      result$pMeans[ i ] <- mean( data[[ pNames[ i ] ]] ) * scalingFactors[ i ]
      result$qMeans[ i ] <- mean( data[[ qNames[ i ] ]] ) / scalingFactors[ i ]
   }
   names( result$pMeans ) <- pNames
   names( result$qMeans ) <- qNames
   if( nFix > 0 ) {
      result$fMeans <- array( NA, nFix )
      for( i in 1:nFix ) {
         result$fMeans[ i ] <- mean( data[[ fNames[ i ] ]] )
      }
      names( result$fMeans ) <- fNames
   }

   ## instrumental variables
   if( nIV == 0 ) {
      inst <- NULL
   } else {
      inst <- as.formula( paste( "~", paste( paste( "iv", c( 1:nIV ), sep="" ),
         collapse = " + " ) ) )
   }

   ## prepare and estimate the model
   modelData <- .snqProfitModelData( data = data, weights = weights,
      pNames = pNames, qNames = qNames, fNames = fNames, ivNames = ivNames,
      form = form, netputScale = scalingFactors, fixedScale = result$fMeans )
   system <- snqProfitSystem( nNetput, nFix )    # equation system
   restrict <- snqProfitRestrict( nNetput, nFix, form )    # restrictions
   result$est <- systemfit( method = method, eqns = system, data = modelData,
      TX = restrict, inst = inst, ... )
   result$coef <- snqProfitCoef( coef = result$est$bt, nNetput = nNetput,
      nFix = nFix, form = form, coefCov = result$est$btcov,
      df = nNetput * nObs - nCoef,
      qNames = qNames, pNames = pNames, fNames = fNames )
      # estimated coefficients
   result$coef <- .snqProfitRescaleCoef( result$coef, nNetput, result$fMeans,
      form )

   result$fitted <- data.frame( profit0 = rep( 0, nObs ) )
   result$residuals <- data.frame( profit0 = rep( 0, nObs ) )
   for( i in 1:nNetput ) {
      result$fitted[[ qNames[ i ] ]] <- result$est$eq[[ i ]]$fitted
      result$fitted[[ "profit0" ]] <- result$fitted[[ "profit0" ]] +
         result$fitted[[ qNames[ i ] ]] * data[[ pNames[ i ] ]] *
         scalingFactors[ i ]
      result$residuals[[ qNames[ i ] ]] <- data[[ qNames[ i ] ]] /
         scalingFactors[ i ] - result$fitted[[ qNames[ i ] ]]
   }
   result$fitted[[ "profit" ]] <- result$fitted[[ "profit0" ]]
   result$fitted[[ "profit0" ]] <- NULL
   result$residuals[[ "profit" ]] <- modelData[[ "profit" ]] -
      result$fitted[[ "profit" ]]
   result$residuals[[ "profit0" ]] <- NULL

   result$r2 <- array( NA, c( nNetput + 1 ) )
   for( i in 1:nNetput ) {
      # result$r2[ i ] <- result$est$eq[[ i ]]$r2
      result$r2[ i ] <- rSquared( data[[ qNames[ i ] ]] / scalingFactors[ i ],
         result$residuals[[ qNames[ i ] ]] )
   }
   result$r2[ nNetput + 1 ] <- rSquared( modelData[[ "profit" ]],
      result$residuals[[ "profit" ]] )
   names( result$r2 ) <- c( qNames, "profit" )

   result$hessian <- snqProfitHessian( result$coef$beta, result$pMeans, weights )
      # Hessian matrix
   result$ela <- snqProfitEla( result$coef$beta, result$pMeans,
      result$qMeans, weights, coefVcov = result$coef$allCoefCov,
      df = result$est$df )   # estimated elasticities
   if( nFix > 0 && form == 0 ) {
      result$fixEla <- snqProfitFixEla( result$coef$delta, result$coef$gamma,
         result$qMeans, result$fMeans, weights )
   }

   result$data     <- data
   result$weights  <- weights
   names( result$weights ) <- pNames
   result$normPrice <- modelData$normPrice
   result$convexity  <- semidefiniteness( result$hessian[
      1:( nNetput - 1 ), 1:( nNetput - 1 ) ] )$positive
   result$pNames  <- pNames
   result$qNames  <- qNames
   result$fNames  <- fNames
   result$ivNames <- ivNames
   result$form    <- form
   result$base    <- base
   result$method  <- method
   result$scalingFactors <- scalingFactors
   names( result$scalingFactors ) <- pNames

   class( result )  <- "snqProfitEst"
   return( result )
}
