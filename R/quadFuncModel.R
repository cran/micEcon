quadFuncModel <- function( yName, xNames, data, shifterNames,
   linear, homWeights, regScale ) {

   checkNames( c( yName, xNames, shifterNames ), names( data ) )

   # check argument 'homWeights'
   .quadFuncCheckHomWeights( homWeights, xNames )

   result <- list()

   result$isPanel <- inherits( data, c( "pdata.frame", "plm.dim" ) )

   if( result$isPanel ) {
      estData <- data[ , 1:2 ]
      estData$y <- data[[ yName ]]
   } else {
      estData <- data.frame( y = data[[  yName ]] )
   }

   if( !is.null( homWeights ) ) {
      estData$deflator <- 0
      for( i in seq( along = homWeights ) ) {
         estData$deflator <- estData$deflator + 
            homWeights[ i ] * data[[ names( homWeights )[ i ] ]]
      }
      xOmit <- names( homWeights )[ 1 ]
      iOmit <- which( xNames == xOmit )
   } else {
      iOmit <- 0
      xOmit <- NULL
   }

   estFormula <- "y ~ 1"
   for( i in seq( along = xNames ) ) {
      if( i != iOmit ) {
         xName <- paste( "a", as.character( i ), sep = "_" )
         estData[[ xName ]] <- .quadFuncVarHom( data, xNames[ i ], 
            homWeights, estData$deflator, xOmit ) / regScale
         estFormula <- paste( estFormula, "+", xName )
      }
   }
   if( !linear ) {
      for( i in seq( along = xNames ) ) {
         for( j in i:length( xNames ) ) {
            if( i != iOmit & j != iOmit ) {
               xName <- paste( "b", as.character( i ), as.character( j ),
                  sep = "_" )
               estData[[ xName ]] <- 0.5 *
                  ifelse( i == j, 1, 2 ) *
                  .quadFuncVarHom( data, xNames[ i ], homWeights, 
                     estData$deflator, xOmit ) * 
                  .quadFuncVarHom( data, xNames[ j ], homWeights, 
                     estData$deflator, xOmit ) / 
                  regScale
               estFormula <- paste( estFormula, "+", xName )
            }
         }
      }
   }
   for( i in seq( along = shifterNames ) ) {
      if( is.factor( data[[ shifterNames[ i ] ]] ) | 
            is.logical( data[[ shifterNames[ i ] ]] ) ) {
         xName <- paste( "d", "_", as.character( i ), "_", sep = "" )
         estData[[ xName ]] <- data[[ shifterNames[ i ] ]]
      } else {
         xName <- paste( "d", as.character( i ), sep = "_" )
         estData[[ xName ]] <- data[[ shifterNames[ i ] ]] / regScale
      }
      estFormula <- paste( estFormula, "+", xName )
   }

   result$estData <- estData
   result$estFormula <- estFormula
   result$iOmit <- iOmit
   
   return( result )
}
