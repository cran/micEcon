quadFuncCalc <- function( xNames, data, coef, shifterNames = NULL,
      homWeights = NULL, quadHalf = TRUE ) {

   # if 'data' is a vector, convert it to a data.frame
   data <- .micEconVectorToDataFrame( data )

   checkNames( c( xNames, shifterNames ), names( data ) )

   # check argument 'homWeights'
   .quadFuncCheckHomWeights( homWeights, xNames )

   nExog <- length( xNames )
   nShifter <- length( shifterNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2 + nShifter

   if( nCoef > length( coef ) ) {
      stop( "a quadratic function with ", nExog, " exogenous variables",
         " and ", nShifter, " shifter variables",
         " must have at least ", nCoef, " coefficients" )
   }

   # calculate index to normalize variables
   if( !is.null( homWeights ) ) {
      deflator <- 0
      for( i in seq( along = homWeights ) ) {
         deflator <- deflator + 
            homWeights[ i ] * data[[ names( homWeights )[ i ] ]]
      }
   }

   result <- rep( coef[ "a_0" ], nrow( data ) )
   for( i in seq( along = xNames ) ) {
      result <- result + coef[ paste( "a", i, sep = "_" ) ] * 
         .quadFuncVarHom( data, xNames[ i ], homWeights, deflator )
      for( j in seq( along = xNames ) ) {
         result <- result + ifelse( quadHalf, 0.5, 1 ) * 
            coef[ paste( "b", min( i, j ), max( i, j ), sep = "_" ) ] *
            .quadFuncVarHom( data, xNames[ i ], homWeights, deflator ) *
            .quadFuncVarHom( data, xNames[ j ], homWeights, deflator )
      }
   }
   for( i in seq( along = shifterNames ) ) {
      if( is.logical( data[[ shifterNames[ i ] ]] ) ) {
         result <- result + coef[ paste( "d", i, "TRUE", sep = "_" ) ] * 
            data[[ shifterNames[ i ] ]]
      } else if( is.factor( data[[ shifterNames[ i ] ]] ) ) {
         for( j in levels( data[[ shifterNames[ i ] ]] ) ) {
            thisCoefName <- paste( "d", i, j, sep = "_" )
            if( thisCoefName %in% names( coef ) ) {
               result <- result + coef[ thisCoefName ] * 
                  ( data[[ shifterNames[ i ] ]] == j )
            }
         }
      } else {
         result <- result + coef[ paste( "d", i, sep = "_" ) ] * 
            data[[ shifterNames[ i ] ]]
      }
   }

   names( result ) <- rownames( data )
   return( result )
}
