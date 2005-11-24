aidsBestA0 <- function( pNames, wNames, xtName,
      data = NULL, ivNames = NULL, method = "IL:L",
      a0min = -50, a0max = 50, stoprange = 3, stopiter = 10,
      verbose = FALSE, ... ) {

   if( length( pNames ) != length( wNames ) ) {
      stop( "arguments 'pNames' and 'wNames' must have the same length" )
   }
   nGoods <- length( pNames )

   if( !( substr( method, 1, 2 ) %in% c( "MK", "IL" ) ) ) {
      stop( "at the moment this function works only with",
         " the 'Iterated Linear Least Squares Estimation' (IL)" )
   }
   if( a0min >= a0max ) stop( "a0min must be smaller than a0max" )

   deta0 <- function( a0, ... ) {
      estResult <- aidsEst( pNames, wNames, xtName, data = data,
         method = method, ivNames = ivNames,
         alpha0 = a0, ... )
      det <- estResult$est$drcov
      assign( "allValues", rbind( allValues, c( a0, det ) ),
         sys.frame( sys.parent( ) )  )
      if( verbose ) {
         cat( "a0:", a0, "-> det:", det, "(iterIL:", estResult$iterIL, ")\n" )
      }
      return( det )
   }

   a0 <- array( NA, c( 2, 4 ) )  # 1st row=alpha0; 2nd row=det(alpha0)
   a0[ 1, ] <- c( a0min, a0min + ( a0max - a0min ) / 3,
      a0max - ( a0max - a0min ) / 3, a0max )
   allValues <- NULL

   for( i in 1:4 ) {
      a0[2,i] <- deta0( a0[ 1, i ], ... )
   }
   iter <- 0
   while( ( which.min( a0[ 2, ] ) == 1 | which.min( a0[ 2, ]) == 4 ) &
      iter < stopiter ) {
      iter <- iter + 1
      if( which.min( a0[ 2, ] ) == 1 ) {
         a0[ , 2:4 ] <- a0[ , 1:3 ]
         a0[ 1, 1 ]  <- a0[ 1, 1 ] - ( a0max - a0min ) / 3
         a0[ 2, 1 ]  <- deta0( a0[ 1, 1 ], ... )
      } else {
         a0[ , 1:3 ] <- a0[ , 2:4 ]
         a0[ 1, 4 ]  <- a0[ 1, 4 ] + ( a0max - a0min ) / 3
         a0[ 2, 4 ]  <- deta0( a0[ 1, 4 ], ... )
      }
   }
   while( iter < stopiter & ( a0[ 1, 4 ] - a0[ 1, 1 ] ) > stoprange ) {
      iter <- iter + 1
      if( which.min( a0[ 2, ] ) == 2 ) {
         a0[ , 4 ] <- a0[ , 3 ]
         if( a0[ 1, 2 ] - a0[ 1, 1 ] >= a0[ 1, 3 ] - a0[ 1, 2 ]) {
            a0[ , 3 ]  <- a0[ , 2 ]
            a0[ 1, 2 ] <- ( a0[ 1, 1 ] + a0[ 1, 3 ] ) / 2
            a0[ 2, 2 ] <- deta0( a0[ 1, 2 ], ... )
         } else {
            a0[ 1, 3 ] <- ( a0[ 1, 2 ] + a0[ 1, 4 ] ) / 2
            a0[ 2, 3 ] <- deta0( a0[ 1, 3 ], ... )
         }
      } else if( which.min(a0[2,])==3) {
         a0[,1] <- a0[,2]
         if( a0[ 1, 4 ] - a0[ 1, 3 ] >= a0[ 1, 3 ] - a0[ 1, 2 ] ) {
            a0[ , 2 ]  <- a0[ , 3 ]
            a0[ 1, 3 ] <- ( a0[ 1, 2 ] + a0[ 1, 4 ] ) / 2
            a0[ 2, 3 ] <- deta0( a0[ 1, 3 ], ... )
         } else {
            a0[ 1, 2 ] <- ( a0[ 1, 1 ] + a0[ 1, 3 ] ) / 2
            a0[ 2, 2 ] <- deta0( a0[ 1, 2 ], ... )
         }
      } else {
         stop("minimum not between a0min and a0max")
      }
   }
   result <- list()
   result$alpha0 <- a0[ 1, which.min( a0[ 2, ] ) ]
   result$allValues <- allValues[ order( allValues[ , 1 ] ), ]
   colnames( result$allValues ) <- c( "a0", "det" )
   result$iter <- iter
   return( result )
}
