micEconIndex <- function( prices, quantities, base, data, method, na.rm,
   na.0, weights, what ) {

   if( length( prices ) != length( quantities ) ) {
      stop( "arguments 'prices' and 'quantities' must have the same length" )
   }

   checkNames( c( prices, quantities ), names( data ) )

   n <- length( prices )

   numerator <- numeric( nrow( data ) )
   denominator <- numeric( nrow( data ) )

   if( method %in% c( "Laspeyres", "Paasche" ) ) {
      for( i in 1:n ) {
         pt <- data[[ prices[ i ] ]]
         p0 <- mean( data[[ prices[ i ] ]][ base ], na.rm = na.rm )
         qt <- data[[ quantities[ i ] ]]
         q0 <- mean( data[[ quantities[ i ] ]][ base ], na.rm = na.rm )
         if( method == "Laspeyres" ) {
            if( is.na( q0 ) || is.na( p0 ) || all( is.na( pt ) ) ) {
               if( !na.0 ) {
                  numerator <- NA
               }
            } else {
               selection <- !is.na( pt )
               numerator[ selection ] <- numerator[ selection ] +
                  pt[ selection ] * q0
               denominator[ selection ] <- denominator[ selection ] + p0 * q0
               if( !na.0 && q0 > 0 ) {
                  numerator[ is.na( pt ) ] <- NA
               }
            }
         } else if( method == "Paasche" ) {
            if( is.na( p0 ) || all( is.na( pt ) ) || all( is.na( qt ) ) ) {
               if( !na.0 ) {
                  numerator <- NA
               }
            } else {
               selection <- qt > 0 & !is.na( qt ) & !is.na( pt )
               numerator[ selection ] <- numerator[ selection ] +
                  pt[ selection ] * qt[ selection ]
               denominator[ selection ] <- denominator[ selection ] +
                  p0 * qt[ selection ]
               if( !na.0 ) {
                  numerator[ is.na( qt ) ] <- NA
                  numerator[ qt > 0  & is.na( pt ) ] <- NA
                  denominator[ is.na( qt ) ] <- NA
               }
            }
         }
      }
      result <- numerator / denominator
      if( weights ) {
         weightData <- data.frame( obsNo = c( 1:nrow( data ) ) )
         rownames( weightData ) <- rownames( data )
         for( i in 1:n ) {
            if( method == "Laspeyres" ) {
               weightData[[ prices[ i ] ]] <- data[[ prices[ i ] ]] *
                  mean( data[[ quantities[ i ] ]][ base ], na.rm = na.rm ) /
                  numerator
            } else if( method == "Paasche" ) {
               weightData[[ prices[ i ] ]] <- data[[ prices[ i ] ]] *
                  data[[ quantities[ i ] ]] / numerator
            }
         }
         weightData$obsNo <- NULL
         attributes( result )$weights <- weightData
      }
   } else if( method == "Fisher" ) {
      pL <- priceIndex( prices, quantities, base, data, method = "Laspeyres",
         na.rm = na.rm, na.0 = na.0, weights )
      pP <- priceIndex( prices, quantities, base, data, method = "Paasche",
         na.rm = na.rm, na.0 = na.0, weights )
      result <- sqrt( pL * pP )
      if( weights ) {
         attributes( result )$weights <-
            0.5 * attributes( pL )$weights +
            0.5 * attributes( pP )$weights
      }
   } else {
      stop( "argument 'method' must be either 'Laspeyres', 'Paasche'",
         " or 'Fisher'" )
   }
   names( result ) <- rownames( data )
   return( result )
}

priceIndex <- function( prices, quantities, base, data,
   method = "Laspeyres", na.rm = FALSE, na.0 = FALSE, weights = FALSE ) {

   checkNames( c( prices, quantities ), names( data ) )

   result <- micEconIndex( prices, quantities, base, data, method, na.rm,
      na.0, weights, "price index" )

   return( result )
}

quantityIndex <- function( prices, quantities, base, data,
   method = "Laspeyres", na.rm = FALSE, na.0 = FALSE, weights = FALSE ) {

   checkNames( c( prices, quantities ), names( data ) )

   result <- micEconIndex( quantities, prices, base, data, method, na.rm,
      na.0, weights, "quantity index" )

   return( result )
}


