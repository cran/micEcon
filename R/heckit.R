heckit <- function( formula, probitformula, data, inst = NULL,
   print.level = 0 ) {

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
   } else if( !is.null( inst ) ) {
      if( class ( inst ) != "formula" || length( inst ) != 2 ) {
         stop( "argument 'inst' must be a 1-sided formula" )
      }
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

   if( print.level > 0 ) {
      cat ( "\nEstimating 1st step Probit model . . ." )
   }
   result$probit <- glm( probitformula, binomial( link = "probit" ), data )
   if( print.level > 0 ) cat( " OK\n" )

   data$probitLambda <- dnorm( result$probit$linear.predictors ) /
      pnorm( result$probit$linear.predictors )

   data$probitDelta <- data$probitLambda * ( data$probitLambda +
      result$probit$linear.predictors )

   step2formula <- as.formula( paste( formula[ 2 ], "~", formula[ 3 ],
      "+ probitLambda" ) )

   if( is.null( inst ) ) {
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step OLS model . . ." )
      }
      result$lm <- lm( step2formula, data, data$probitdummy == 1 )
      resid <- residuals( result$lm )
       step2coef <- coefficients( result$lm )
      if( print.level > 0 ) cat( " OK\n" )
   } else {
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step 2SLS/IV model . . ." )
      }
      formulaList <- list( step2formula )
      instImr <- as.formula( paste( "~", inst[ 2 ], "+ probitLambda" ) )
      library( systemfit )
      result$lm <- systemfit( "2SLS", formulaList, inst = instImr,
         data = data[ data$probitdummy == 1, ] )
      resid <- residuals( result$lm )[ , 1 ]
       step2coef <- coefficients( result$lm$eq[[ 1 ]] )
      if( print.level > 0 ) cat( " OK\n" )
   }

   result$sigma <- as.numeric( sqrt( crossprod( resid ) /
      sum( data$probitdummy == 1 ) +
      mean( data$probitDelta[ data$probitdummy == 1 ] ) *
       step2coef[ "probitLambda" ]^2 ) )

   result$rho <-  step2coef[ "probitLambda" ] / result$sigma
   result$probitLambda <- data$probitLambda
   result$probitDelta  <- data$probitDelta
   if( print.level > 0 ) {
      cat ( "Calculating coefficient covariance matrix . . ." )
   }
   # the foolowing variables are named according to Greene (2003), p. 785
   if( is.null( inst ) ) {
      xMat <- model.matrix( result$lm )
   } else {
      xMat <- result$lm$eq[[ 1 ]]$x
   }
   wMat <- model.matrix( result$probit )[ data$probitdummy == 1, ]
   #fMat <- t( xMat ) %*% diag( result$probitDelta[
   #   data$probitdummy == 1 ] ) %*% wMat
   # replaced the previous lines by the following to avoid the
   # diagonal matrix that gets too large for large data sets.
   txdMat <- t( xMat )
   dVec <- result$probitDelta[  data$probitdummy == 1 ]
   for( i in 1:nrow( txdMat ) ) {
      txdMat[ i, ] <- txdMat[ i, ] * dVec
   }
   fMat <- txdMat %*% wMat
   rm( txdMat, dVec )
   qMat <- result$rho^2 * ( fMat %*% vcov( result$probit )%*% t( fMat ) )
   #result$vcov<-result$sigma^2*solve(crossprod(xMat))%*%
   #(t(xMat)%*%diag(1-result$rho^2*
   #result$probitDelta[data$probitdummy==1])%*%
   #(txd2Mat%*%
   #xMat+qMat)%*%solve(crossprod(xMat))
   # replaced the previous lines by the following to avoid the
   # diagonal matrix that gets too large for large data sets.
   txd2Mat <- t( xMat )
   d2Vec <-  1 - result$rho^2 * result$probitDelta[ data$probitdummy == 1 ]
   for( i in 1:nrow( txd2Mat ) ) {
      txd2Mat[ i, ] <- txd2Mat[ i, ] * d2Vec
   }
   result$vcov <- result$sigma^2 * solve( crossprod( xMat ) ) %*%
      ( txd2Mat %*%
      xMat + qMat ) %*% solve( crossprod( xMat ) )
   rm( txd2Mat, d2Vec )
   if( print.level > 0 ) cat( " OK\n" )
   result$coef <- matrix( NA, nrow = length( step2coef ), ncol = 4 )
   rownames( result$coef ) <- names( step2coef )
   colnames( result$coef ) <- c( "Estimate", "Std. Error", "t value",
      "Pr(>|t|)" )
   result$coef[ , 1 ] <- step2coef
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

   if( class( x$lm ) == "lm" ) {
      rSquared <- c( summary( x$lm )$r.squared, summary( x$lm )$adj.r.squared )
   } else {
      rSquared <- c( x$lm$eq[[ 1 ]]$r2, x$lm$eq[[ 1 ]]$adjr2 )
   }
   cat( paste(
      "Multiple R-Squared:", round( rSquared[ 1 ], digits),
      "Adjusted R-Squared:", round( rSquared[ 2 ], digits),
      "\n" ) )
   cat("\n")
   invisible( x )
}
