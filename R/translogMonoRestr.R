translogMonoRestr <- function( xNames, data, quadHalf = TRUE,
   logValues = FALSE, area = FALSE ) {

   checkNames( c( xNames ), names( data ) )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( logValues ) {
      logData <- data
   } else {
      logData   <- data.frame( no = c( 1:nrow( data ) ) )
      for( i in seq( along = xNames ) ) {
         logData[[ xNames[ i ] ]] <- log( data[[ xNames[ i ] ]] )
      }
   }

   if( area ) {
      extremeLogValues <- list()
      for( i in seq( along = xNames ) ) {
         extremeLogValues[[ i ]] <- c(
            min( logData[[ xNames[ i ] ]] ),
            max( logData[[ xNames[ i ] ]] ) )
      }
      logData <- expand.grid( extremeLogValues )
      colnames( logData ) <- xNames
   }
   nObs  <- nrow( logData )
   restr <- matrix( 0, nObs * nExog, nCoef )
   for( i in seq( along = xNames ) ) {
      myRows <- c( ( ( i - 1 ) * nObs + 1 ):( i * nObs ) )
      restr[ myRows, 1 + i ] <- 1
      for( j in seq( along = xNames ) ) {
         restr[ myRows, 1 + nExog + veclipos( i, j, nExog ) ] <-
            ifelse( quadHalf, 1, 2 ) * logData[[ xNames[ j ] ]]
      }
   }

   return( restr )
}

# test with (only if area == FALSE):
# matrix( translogMonoRestr( lnInputNames, estData ) %*% a$allCoef, ncol=4) *
# exp(translogCalc( lnInputNames, estData, a$allCoef) %*% t(c(1,1,1,1))) /
# estData[ , inputNames ] -
# translogDeriv( lnInputNames, estData, a$allCoef )$deriv