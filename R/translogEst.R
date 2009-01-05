translogEst <- function( yName, xNames, data, shifterNames = NULL,
   quadHalf = TRUE, dataLogged = FALSE, ... ) {

   checkNames( c( yName, xNames, shifterNames ), names( data ) )

   if( dataLogged ) {
      logData   <- data
   } else {
      if( "plm.dim" %in% class( data ) ) {
         logData <- data[ , 1:2 ]
      } else {
         logData <- data.frame( no = c( 1:nrow( data ) ) )
      }
      logData[[ yName ]] <- log( data[[ yName ]] )
      for( i in seq( along = xNames ) ) {
         logData[[ xNames[ i ] ]] <- log( data[[ xNames[ i ] ]] )
      }
      for( i in seq( along = shifterNames ) ) {
         logData[[ shifterNames[ i ] ]] <-
            log( data[[ shifterNames[ i ] ]] )
      }
   }

   result <- quadFuncEst( yName = yName, xNames = xNames, data = logData,
      shifterNames = shifterNames, quadHalf = quadHalf, ... )

   result$r2nonLog <- rSquared( exp( logData[[ yName ]] ),
      exp( logData[[ yName ]] ) - exp( result$fitted ) )

   if( !dataLogged ){
      result$fitted <- exp( result$fitted )
   }

   class( result ) <- "translogEst"
   return( result )
}
