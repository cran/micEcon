heckit <- function( formula, probitformula, data ) {

   if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   } else if( length( formula ) != 3 ) {
      stop( "argument 'formula' must be a 2-sided formula" )
   } else if( "probit" %in% substr( all.vars( formula ), 1, 6 ) ) {
      stop( paste( "argument 'formula' may not include variable names",
      "starting with 'probit'" ) )
   } else if( class( probitformula ) != "formula" ) {
      stop( "argument 'probitformula' must be a formula" )
   } else if( length( probitformula ) != 3 ) {
      stop( "argument 'probitformula' must be a 2-sided formula" )
   } else if( "probit" %in% substr( all.vars( probitformula ), 1, 6 ) ) {
      stop( paste( "argument 'probitformula' may not include a variable",
         "names starting with 'probit'" ) )
   }

   result <- list()

   data$probitdummy <- model.frame( probitformula, data = data )[ , 1 ]
   test <- levels( as.factor( as.numeric( data$probitdummy ) ) )
   if( length( test ) != 2 ) {
      stop( paste( "The left hand side of 'probitformula' may only contain",
         "1 and 0 or TRUE and FALSE" ) )
   } else if( !all.equal( test, c( "0", "1" ) ) ) {
      stop( paste( "The left hand side of 'probitformula' may only contain",
         "1 and 0 or TRUE and FALSE" ) )
   }

   result$probit <- glm( probitformula, binomial( link = "probit" ), data )

   data$probitLambda <- dnorm( result$probit$linear.predictors ) /
      pnorm( result$probit$linear.predictors )

   data$probitDelta <- data$probitLambda * ( data$probitLambda +
      result$probit$linear.predictors )

   step2formula <- as.formula( paste( formula[ 2 ], "~", formula[ 3 ],
      "+ probitLambda" ) )

   result$lm <- lm( step2formula, data, data$probitdummy == 1 )

   result$sigma <- as.numeric( sqrt( crossprod( residuals( result$lm ) ) /
      sum( data$probitdummy == 1 ) +
      mean( data$probitDelta[ data$probitdummy == 1 ] ) *
      coefficients( result$lm )[ "probitLambda" ]^2 ) )

   result$rho <- coefficients( result$lm )[ "probitLambda" ] / result$sigma
   result$probitLambda <- data$probitLambda
   result$probitDelta  <- data$probitDelta

   # the foolowing variables are named according to Greene (2003), p. 785
   xMat <- model.matrix( result$lm )
   wMat <- model.matrix( result$probit )[ data$probitdummy == 1, ]
   fMat <- t( xMat ) %*% diag( result$probitDelta[
      data$probitdummy == 1 ] ) %*% wMat
   qMat <- result$rho^2 * ( fMat %*% vcov( result$probit )%*% t( fMat ) )
   result$vcov <- result$sigma^2 * solve( crossprod( xMat ) ) %*%
      ( t( xMat ) %*% diag( 1 - result$rho^2 *
      result$probitDelta[ data$probitdummy == 1 ] ) %*%
      xMat + qMat ) %*% solve( crossprod( xMat ) )
   result$coef <- matrix( NA, nrow = length( coef( result$lm ) ), ncol = 4 )
   rownames( result$coef ) <- names( coef( result$lm ) )
   colnames( result$coef ) <- c( "Estimate", "Std. Error", "t value",
      "Pr(>|t|)" )
   result$coef[ , 1 ] <- coef( result$lm )
   result$coef[ , 2 ] <- sqrt( diag( result$vcov ) )
   result$coef[ , 3 ] <- result$coef[ , 1 ] / result$coef[ , 2 ]
   result$coef[ , 4 ] <- 2 * ( 1 - pt( abs( result$coef[ , 3 ] ),
      result$lm$df ) )

   class( result ) <- "heckit"
   return( result )
}

summary.heckit <- function( object, ... ) {
   print( object, ... )
   invisible( object )
}

print.heckit <- function( x, digits = 6, ... ) {
   Signif <- symnum( x$coef[ , 4 ], corr = FALSE, na = FALSE,
      cutpoints = c( 0, 0.001, 0.01, 0.05, 0.1, 1 ),
      symbols   = c( "***", "**", "*", "." ," " ))

   table <- cbind( round( x$coef, digits ), Signif )

   rownames( table ) <- rownames( x$coef )
   colnames( table ) <- c( "Estimate", "Std. Error", "t value",
      "Pr(>|t|)", "" )

   print.matrix( table, quote = FALSE, right = TRUE )
   cat( "---\nSignif. codes: ", attr( Signif, "legend" ), "\n" )

   cat( paste(
      "Multiple R-Squared:", round( summary( x$lm )$r.squared, digits),
      "Adjusted R-Squared:", round( summary( x$lm )$adj.r.squared, digits),
      "\n" ) )
   cat("\n")
   invisible( x )
}
