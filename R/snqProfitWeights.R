## ---- snqProfit: default weights for normalizing prices ----------
snqProfitWeights <- function( pNames, qNames, data, method = "DW92", base = 1 ) {

   if( length( qNames ) != length( pNames ) ) {
      stop( "arguments 'qNames' and 'pNames' must have the same length." )
   }

   nNetput <- length( qNames )
   weights <- array( 0, c( nNetput ) )   # weights of the netput prices

   if( method == "DW92" ) {
      totalValue <- 0
      for( i in 1:nNetput ) {
         totalValue <- totalValue + mean( abs( data[[ qNames[ i ] ]] ) ) *
            mean( data[[ pNames[ i ] ]][ base ] )
      }
      for( i in 1:nNetput ) {
         weights[ i ] <- mean( abs( data[[ qNames[ i ] ]] ) ) *
            mean( data[[ pNames[ i ] ]][ base ] ) / totalValue
      }
   } else if( method == "Kohli" ) {
      totalValue <- 0
      for( i in 1:nNetput ) {
         totalValue <- totalValue + mean( data[[ pNames[ i ] ]] ) *
            mean( abs( data[[ qNames[ i ] ]] ) )
      }
      for( i in 1:nNetput ) {
         weights[ i ] <-
            mean( abs( data[[ qNames[ i ] ]] ) ) /
            totalValue
      }
   } else if( method == "Ooms" ) {
      totalValues <- 0
      for( i in 1:nNetput ) {
         totalValues <- totalValues +
            abs( data[[ pNames[ i ] ]] * data[[ qNames[ i ] ]] )
      }
      for( i in 1:nNetput ) {
         weights[ i ] <- mean( abs( data[[ pNames[ i ] ]] *
            data[[ qNames[ i ] ]] ) ) / mean( totalValues )
      }
   } else {
      stop( "argument 'message' must be either 'DW92', 'Kohli' or 'Ooms'." )
   }

   return( weights )
}
