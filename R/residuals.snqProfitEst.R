residuals.snqProfitEst <- function( object, scaled = TRUE, ... ) {

   nNetput <- length( object$pMeans )
   nFixed  <- length( object$fMeans )
   nObs    <- nrow( object$data )

   result <- data.frame( profit0 = rep( 0, nObs ) )
   for( i in 1:nNetput ) {
      result[[ object$qNames[ i ] ]] <-
         object$data[[ object$qNames[ i ] ]] /
            object$scalingFactors[ i ]^( scaled ) -
         object$fitted[[ object$qNames[ i ] ]] *
            object$scalingFactors[ i ]^( !scaled )
      result$profit0 <- result$profit0 +
         object$data[[ object$qNames[ i ] ]] *
         object$data[[ object$pNames[ i ] ]]
   }
   result$profit <- result$profit0 - object$fitted$profit
   result$profit0 <- NULL

   return( result )
}

## the same for snqProfitImposeConvexity
residuals.snqProfitImposeConvexity <- function( object, scaled = TRUE, ... ) {

   result <- residuals.snqProfitEst( object, scaled = scaled, ... )

   return( result )
}
