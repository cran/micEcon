translogProdFuncMargCost <- function( yName, xNames, wNames,
      data, coef, dataLogged = FALSE ) {

   checkNames( c( yName, xNames, wNames ), names( data ) )

   if( length( yName ) != 1 ) {
      stop( "argument 'yName' must include the name of one single output" )
   }

   if( length( xNames ) != length( wNames ) ) {
      stop( "arguments 'xNames' and 'wNames' must have the same length" )
   }

   if( dataLogged ) {
      logData   <- data
   } else {
      logData <- logDataSet( data = data,
         varNames = c( yName, xNames, wNames ) )
   }

   # number of inputs
   nInput <- length( xNames )

   # coefficients with names as if theta was a normal input
   coefTl <- coef

   # calculate distance and theta
   xNamesTheta <- xNames

   # elasticities: d log F / d log x_i
   elaData <- translogEla( xNames = xNamesTheta, data = logData, coef = coefTl,
      dataLogged = TRUE )

   # matrix of beta coefficients
   beta <- matrix( NA, nrow = nInput, ncol = nInput )
   for( i in 1:nInput ) {
      for( j in i:nInput ) {
         beta[ i, j ] <- coef[ paste( "b", i, j, sep = "_" ) ]
         beta[ j, i ] <- beta[ i, j ]
      }
   }

   # marginal costs
   margCost <- rep( NA, nrow( logData ) )

   for( t in 1:nrow( logData ) ) {
      yVal <- exp( as.numeric( logData[ t, yName ] ) )
      xVal <- exp( as.numeric( logData[ t, xNames ] ) )
      wVal <- exp( as.numeric( logData[ t, wNames ] ) )
      eVal <- as.numeric( elaData[ t, ] )

      # Jacobian matrix of g with respect to x
      gxJac <- matrix( NA, nrow = nInput, ncol = nInput )
      for( j in 1:nInput ) {
         for( i in 1:( nInput - 1 ) ) {
            gxJac[ i, j ] <-
               ( i == j ) * wVal[ nInput ] *
                  ( xVal[ nInput ] / xVal[ j ]^2 ) *
                  eVal[ i ] / eVal[ nInput ] -
               ( j == nInput ) * wVal[ nInput ] * ( 1 / xVal[ i ] ) *
                  eVal[ i ] / eVal[ nInput ] -
               wVal[ nInput ] * ( xVal[ nInput ] / xVal[ i ] ) *
                  ( beta[ i, j ] / xVal[ j ] ) / eVal[ nInput ] +
               wVal[ nInput ] * ( xVal[ nInput ] / xVal[ i ] ) *
                  ( eVal[ i ] / eVal[ nInput ]^2 ) *
                  beta[ nInput, j ] / xVal[ j ]
         }
         gxJac[ nInput, j ] <- eVal[ j ] / xVal[ j ]
      }

      # Jacobian matrix of g with respect to y
      gyJac <- rep( NA, nInput )
      gyJac[ 1:( nInput - 1 ) ] <- 0
      gyJac[ nInput ] <- - 1 / yVal

      # Jacobian matrix of x with respect to y
      xyJac <- - solve( gxJac, gyJac )

      # marginal costs
      margCost[ t ] <- t( wVal ) %*% xyJac
   }

   names( margCost ) <- rownames( data )

   return( margCost )
}