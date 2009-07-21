cesEst <- function( yName, xNames, data, vrs = FALSE, ... ) {

   # y = gamma * ( alpha * x1^rho + ( 1 - alpha ) * x2^rho )^(phi/rho)
   # s = 1 / ( 1 - rho )

   checkNames( c( yName, xNames ), names( data ) )

   if( length( xNames ) != 2 ) {
      stop( "currently, argument 'xNames' must contain exactly",
         " two variable names" )
   }

   # store the (matched) call
   matchedCall <- match.call()

   # prepare data for estimation
   estData <- data.frame( y = data[[ yName ]],
      x1 = data[[ xNames[ 1 ] ]], x2 = data[[ xNames[ 2 ] ]] )

   # start values
   startVal <- c( gamma = sqrt( mean( estData$y ) ),
      alpha = 0.5, rho = 0.5 )
   if( vrs ) {
      startVal <- c( startVal, phi = 1 )
   }

   cesRss <- function( par, data ) {
      gamma <- par[ "gamma" ]
      alpha <- par[ "alpha" ]
      rho <- par[ "rho" ]
      if( "phi" %in% names( par ) ) {
         phi <- par[ "phi" ]
      } else {
         phi <- 1
      }
      yHat <- gamma *
         ( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho )^( phi / rho )
      return( sum( ( data$y - yHat )^2 ) )
   }

   cesRssDeriv <- function( par, data ) {
      result <- par
      gamma <- par[ "gamma" ]
      alpha <- par[ "alpha" ]
      rho <- par[ "rho" ]
      if( "phi" %in% names( par ) ) {
         phi <- par[ "phi" ]
      } else {
         phi <- 1
      }
      yHat <- gamma *
         ( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho )^( phi / rho )
      resid <- data$y - yHat
      result[ "gamma" ] <- sum( - 2 * resid * 
         ( ( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho )^( phi / rho ) ) )
      result[ "alpha" ] <- sum( -2 * resid * ( gamma * ( phi / rho ) *
         ( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho )^( phi / rho - 1 ) *
         ( data$x1^rho - data$x2^rho ) ) )
      result[ "rho" ] <- sum( -2 * resid * ( gamma *
         log( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho ) *
         ( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho )^( phi / rho ) *
         ( - phi / rho^2 ) +
         gamma * ( phi / rho ) *
         ( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho )^( phi / rho - 1 ) *
         ( alpha * log( data$x1 ) * data$x1^rho +
            ( 1 - alpha ) * log( data$x2 ) * data$x2^rho ) ) )
      if( "phi" %in% names( par ) ) {
         result[ "phi" ] <- sum( -2 * resid * ( gamma *
            log( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho ) *
            ( alpha * data$x1^rho + ( 1 - alpha ) * data$x2^rho )^( phi / rho ) *
            ( 1 / rho ) ) )
      }

      return( result )
   }

   if( is.null( matchedCall$method ) || matchedCall$method == "SANN" ) {
      result <- optim( par = startVal, fn = cesRss, data = estData,
         hessian = TRUE, ... )
   } else {
      result <- optim( par = startVal, fn = cesRss, gr = cesRssDeriv,
         data = estData, hessian = TRUE, ... )
   }

   # covariance matrix of the estimated parameters
   if( det( result$hessian ) >= .Machine$double.eps ) {
      result$vcov <- solve( result$hessian )
   } else {
      result$vcov <- matrix( NA, nrow = length( result$par ),
         ncol = length( result$par ) )
      dimnames( result$vcov ) <- dimnames( result$hessian )
   }

   # return also the call
   result$call <- matchedCall 

   # nonlinear least squares
#    result$nls <- nls(
#       y ~ gamma * ( alpha * x1^rho + ( 1 - alpha ) * x2^rho )^(1/rho),
#       data = estData, start = result$startVal, trace = TRUE,
#       algorithm = "port", lower = c( -Inf, 0.01 , -Inf ),
#       upper = c( Inf, 0.99, Inf ) )

   class( result ) <- c( "cesEst", class( result ) )
   return( result )
}

