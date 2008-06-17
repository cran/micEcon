quadFuncEst <- function( yName, xNames, data, quadHalf = TRUE, exVarScale = 1 ) {

   checkNames( c( yName, xNames ), names( data ) )

   nExog   <- length( xNames )
   result <- list()

   estData <- data.frame( y = data[[  yName ]] )
   estFormula <- "y ~ 1"
   for( i in 1:nExog ) {
      xName <- paste( "x", as.character( i ), sep = "" )
      estData[[ xName ]] <- data[[ xNames[ i ] ]] / exVarScale
      estFormula <- paste( estFormula, "+", xName )
   }
   for( i in 1:nExog ) {
      for( j in i:nExog ) {
         xName <- paste( "x", as.character( i ), ".", as.character( j ),
            sep = "" )
         estData[[ xName ]] <- ifelse( quadHalf, 0.5, 1 ) *
            ifelse( i == j, 1, 2 ) *
            data[[ xNames[ i ] ]] * data[[ xNames[ j ] ]] / exVarScale
         estFormula <- paste( estFormula, "+", xName )
      }
   }
   result$lm <- lm( as.formula( estFormula ), estData )
   result$residuals <- result$lm$residuals
   result$fitted    <- result$lm$fitted.values

   # coefficients and their covariance matrix
   result$coef      <- coef( result$lm )
   result$coefCov   <- vcov( result$lm )
   coefNames <- .quadFuncCoefNames( nExog )
   names( result$coef )      <- coefNames
   rownames( result$coefCov ) <- coefNames
   colnames( result$coefCov ) <- coefNames

   result$r2    <- summary( result$lm )$r.squared
   result$r2bar <- summary( result$lm )$adj.r.squared
   result$nObs  <- length( result$lm$residuals )
   result$model.matrix <- cbind( rep( 1, result$nObs ),
      as.matrix( estData[ , 2:( ncol( estData ) ) ] ) )
   class( result ) <- "quadFuncEst"
   return( result )
}
