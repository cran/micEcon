heckit <- function( selection, formula, data, inst = NULL,
   na.action = options( "na.action" ), print.level = 0 ) {

   if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   } else if( length( formula ) != 3 ) {
      stop( "argument 'formula' must be a 2-sided formula" )
   } else if( "invMillsRatio" %in% all.vars( formula ) ) {
      stop( "argument 'formula' may not include a variable name",
         " 'invMillsRatio'" )
   } else if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   } else if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   } else if( "probit" %in% substr( all.vars( selection ), 1, 6 ) ) {
      stop( "argument 'selection' may not include a variable",
         " names starting with 'probit'" )
   } else if( !is.null( inst ) ) {
      if( class ( inst ) != "formula" || length( inst ) != 2 ) {
         stop( "argument 'inst' must be a 1-sided formula" )
      }
   }

   result <- list()

   probitEndogenous <- model.frame( selection, data = data,
      na.action = NULL )[ , 1 ]
   probitLevels <- levels( as.factor( probitEndogenous ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   probitDummy <- probitEndogenous == probitLevels[ 2 ]

   # NA action
   firstStepData <- model.frame( selection, data = data, na.action = NULL )
   secondStepData <- model.frame( formula, data = data, na.action = NULL )
   if( !is.null( inst ) ) {
      secondStepData <- cbind( secondStepData,
         model.frame( inst, data = data, na.action = NULL ) )
   }
   firstStepOk <- rowSums( is.na( firstStepData ) ) == 0
   secondStepOk <- rowSums( is.na( secondStepData ) ) == 0
   result$dataOk <- firstStepOk & ( secondStepOk | !probitDummy )

   if( na.action == "na.omit" ) {
      data <- data[ result$dataOk, ]
      probitDummy <- probitDummy[ result$dataOk ]
   } else if( na.action == "na.fail" ) {
      if( sum( !firstStepOk ) > 0 ) {
         stop( "missing values at the first step" )
      }
      if( sum( !( secondStepOk | !probitDummy ) ) > 0 ) {
         stop( "missing values at the second step" )
      }
   }

   if( print.level > 0 ) {
      cat ( "\nEstimating 1st step Probit model . . ." )
   }
#   result$probit <- glm( selection, binomial( link = "probit" ), data )
   result$probit <- probit(selection, data=data, x=TRUE, print.level=print.level - 1)
   if( print.level > 0 ) {
       cat( " OK\n" )
   }
#    data$probitLambda <- dnorm(linearPredictors(result$probit)) /
#        pnorm(linearPredictors(result$probit))
#    data$probitDelta <- data$probitLambda * ( data$probitLambda +
#                                             linearPredictors(result$probit))
   imrData <- invMillsRatio( result$probit )
   data$invMillsRatio <- imrData$IMR1
   result$imrDelta <- imrData$delta1

   step2formula <- as.formula( paste( formula[ 2 ], "~", formula[ 3 ],
                                     "+ invMillsRatio" ) )

   if( is.null( inst ) ) {
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step OLS model . . ." )
      }
      result$lm <- lm( step2formula, data, subset = probitDummy )
      resid <- residuals( result$lm )
       step2coef <- coefficients( result$lm )
      if( print.level > 0 ) cat( " OK\n" )
   } else {
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step 2SLS/IV model . . ." )
      }
      formulaList <- list( step2formula )
      instImr <- as.formula( paste( "~", inst[ 2 ], "+ invMillsRatio" ) )
      library( systemfit )
      result$lm <- systemfit( "2SLS", formulaList, inst = instImr,
         data = data[ probitDummy, ] )
      resid <- residuals( result$lm )[ , 1 ]
       step2coef <- coefficients( result$lm$eq[[ 1 ]] )
      if( print.level > 0 ) cat( " OK\n" )
   }
   result$sigma <- as.numeric( sqrt( crossprod( resid ) /
      sum( probitDummy ) +
      mean(result$imrDelta[ probitDummy ] ) *
       step2coef[ "invMillsRatio" ]^2 ) )
   result$rho <-  step2coef[ "invMillsRatio" ] / result$sigma
   names(result$rho) <- NULL
                                        # otherwise the name of step2coef is left...
   result$invMillsRatio <- data$invMillsRatio
   if( print.level > 0 ) {
      cat ( "Calculating coefficient covariance matrix . . ." )
   }
   # the following variables are named according to Greene (2003), p. 785
   if( is.null( inst ) ) {
      xMat <- model.matrix( result$lm )
   } else {
      xMat <- result$lm$eq[[ 1 ]]$x
   }
   result$vcov <- heckitVcov( xMat,
      model.matrix( result$probit )[ probitDummy, ],
      vcov( result$probit ),
      result$rho,
      result$imrDelta[ probitDummy ],
      result$sigma )
# result$vcov <- result$sigma^2 * solve( crossprod( xMat ) ) %*%
#      ( txd2Mat %*%
#      xMat + qMat ) %*% solve( crossprod( xMat ) )
#   rm( txd2Mat, d2Vec )
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

   print( table, quote = FALSE, right = TRUE )
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
