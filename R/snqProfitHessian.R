## ------- calculation of Hessian -------------
snqProfitHessian <- function( beta, prices, weights ) {
   prices <- unlist( prices )
   normPrice <- sum( t( prices ) %*% weights )
   Hessian <- beta / normPrice -
      beta %*% prices %*% t( weights ) / normPrice^2 -
      weights %*% t( prices ) %*% beta / normPrice^2 +
      weights %*% t( weights ) *
      mean( ( t( prices ) %*% beta %*% prices ) / normPrice^3 )
   return( Hessian )
}
