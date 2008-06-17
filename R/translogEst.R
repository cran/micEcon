translogEst <- function( yName, xNames, data, quadHalf = TRUE,
   dataLogged = FALSE ) {

   checkNames( c( yName, xNames ), names( data ) )

   if( dataLogged ) {
      logData   <- data
   } else {
      logData <- data.frame( no = c( 1:nrow( data ) ) )
      logData[[ yName ]] <- log( data[[ yName ]] )
      for( i in seq( along = xNames ) ) {
         logData[[ xNames[ i ] ]] <- log( data[[ xNames[ i ] ]] )
      }
   }

   result <- quadFuncEst( yName, xNames, logData, quadHalf = quadHalf )

   result$r2nonLog <- rSquared( exp( logData[[ yName ]] ),
      exp( logData[[ yName ]] ) - exp( result$fitted ) )

   if( !dataLogged ){
      result$fitted <- exp( result$fitted )
   }

   class( result ) <- "translogEst"
   return( result )
}
