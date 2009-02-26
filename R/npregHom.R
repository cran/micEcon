npregHom <- function( yName, xNames, homWeights, data,
   restrictGrad = TRUE, bws = NULL, ... ) {

   checkNames( c( yName, xNames ), names( data ) )

   # check argument 'homWeights'
   .quadFuncCheckHomWeights( homWeights, xNames )

   result <- list()
   result$call <- match.call()
   result$yName  <- yName
   result$xNames <- xNames
   result$homWeights <- homWeights
   result$restrictGrad <- restrictGrad

   # endogenous variable
   yData <- data[[ yName ]]

   # variable for normalization
   weightVarMatrix <- matrix( NA, nrow = nrow( data ),
      ncol = length( homWeights ) )
   for( i in seq( along = homWeights ) ) {
      varName <- names( homWeights )[ i ]
      if( "factor" %in% class( data[[ varName ]] ) ) {
         stop( "homogeneity can be imposed only on continuous variables",
            " but variable '", varName, "' is an (un)ordered factor" )
      }
      weightVarMatrix[ , i ] <- homWeights[ i ] * data[[ varName ]]
   }
   deflator <- rowSums( weightVarMatrix )
   whichHom <- which( xNames %in% names( homWeights ) )
   if( restrictGrad ) {
      xOmit <- names( homWeights )[ 1 ]
      iOmit <- which( xNames == xOmit )
   } else {
      xOmit <- NULL
      iOmit <- 0
   }

   # regressors
   xData <- data.frame( no = 1:nrow( data ) )
   for( i in seq( along = xNames ) ) {
      if( i != iOmit ) {
         xName <- paste( "r", as.character( i ), sep = "_" )
         xData[[ xName ]] <- .quadFuncVarHom( data, xNames[ i ],
            homWeights, deflator, xSubtract = xOmit )
      }
   }
   xData$no <- NULL

   # nonparametric regression
   if( is.null( bws ) ) {
      result$est <- npreg( tydat = yData, txdat = xData, gradients = TRUE, ... )
   } else {
      result$est <- npreg( tydat = yData, txdat = xData, gradients = TRUE,
         bws = bws, ... )
   }

   # gradients
   if( restrictGrad ) {
      npAllGrad <- cbind( rep( 0, nrow( result$est$grad ) ), result$est$grad )
      colnames( npAllGrad ) <- c( paste( "r", iOmit, sep = "_" ),
         colnames( xData ) )
   } else {
      npAllGrad <- result$est$grad
      colnames( npAllGrad ) <- colnames( xData )
   }
   for( i in whichHom[ whichHom != iOmit ] ) {
      npAllGrad[ , 1 ] <- npAllGrad[ , 1 ] -
         npAllGrad[ , paste( "r", i, sep = "_" ) ]
   }
   result$grad <- matrix( NA, nrow = nrow( xData ), ncol = ncol( npAllGrad ) )
   colnames( result$grad ) <- xNames
   for( i in seq( along = xNames ) ) {
      result$grad[ , i ] <- npAllGrad[ , paste( "r", i, sep = "_" ) ]
      if( i %in% whichHom ) {
         result$grad[ , i ] <- result$grad[ , i ] / deflator
         for( j in whichHom ) {
            result$grad[ , i ] <- result$grad[ , i ] -
               homWeights[ xNames[ i ] ] *
               npAllGrad[ , paste( "r", j, sep = "_" ) ] *
               data[[ xNames[ j ] ]] / deflator^2
         }
      }
   }

   class( result ) <- "npregHom"
   return( result )
}