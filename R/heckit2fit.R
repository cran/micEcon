heckit2fit <- function( selection, outcome,
                   data=sys.frame(sys.parent()),
                   inst = NULL,
                   print.level = 0, ... ) {
   # What is the role of na.action here?  We cannot use na.omit -- we must not omit the observation
   # where outcome is not observed.  na-s cannot be passed either.
   # However, we can (and should?) omit the na-s in explanatory and probit outcomes.  This needs
   # a bit of refinement.
   if( class( outcome ) != "formula" ) {
      stop( "argument 'outcome' must be a formula" )
   } else if( length( outcome ) != 3 ) {
      stop( "argument 'outcome' must be a 2-sided formula" )
   } else if( "invMillsRatio" %in% all.vars( outcome ) ) {
      stop( "argument 'outcome' may not include a variable name",
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
   result$call <- match.call()
   ## Now extract model frames etc.
   ## Selection equation
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset", "weights", 
                "offset"), names(mf), 0)
   mfS <- mf[c(1, m)]
   mfS$na.action <- na.pass
   mfS$drop.unused.levels <- TRUE
   mfS[[1]] <- as.name("model.frame")
   names(mfS)[2] <- "formula"
                                        # model.frame requires the parameter to
                                        # be 'formula'
   mfS <- eval(mfS, parent.frame())
   mtS <- attr(mfS, "terms")
   XS <- model.matrix(mtS, mfS)
   NXS <- ncol(XS)
   YS <- factor(model.response(mfS, "numeric"))
   badRow <- is.na(YS)
   badRow <- badRow | apply(XS, 1, function(v) any(is.na(v)))
                                        # check for NA-s.  Because we have to find NA-s in several
                                        # frames, we cannot use the standard na.* functions here.
                                        # Find bad rows and remove them later.
   probitLevels <- levels( as.factor( YS ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   probitDummy <- YS == probitLevels[ 2 ]
   ## Outcome equation
   m <- match(c("outcome", "data", "subset", "weights",
                "offset"), names(mf), 0)
   mfO <- mf[c(1, m)]
   mfO$na.action <- na.pass
   mfO$drop.unused.levels <- TRUE
   mfO$na.action <- na.pass
                                        # Here we have to keep NA-s: unobserved outcome variables may
                                        # be marked as NA-s.
   mfO[[1]] <- as.name("model.frame")
                                        # eval it as model frame
   names(mfO)[2] <- "formula"
   mfO <- eval(mfO, parent.frame())
   mtO <- attr(mfO, "terms")
   XO <- model.matrix(mtO, mfO)
                                        # the explanatory variables in matrix form
   NXO <- ncol(XO)
   YO <- model.response(mfO, "numeric")
   ## Remove NA observations
   badRow <- badRow | (is.na(YO) & (!is.na(YS) & YS == 1))
   badRow <- badRow | (apply(XO, 1, function(v) any(is.na(v))) & (!is.na(YS) & YS == 1))
                                        # rows in outcome, which contain NA and are observable -> bad too
   if(print.level > 0) {
      cat(sum(badRow), "invalid observations\n")
   }
   XS <- XS[!badRow,]
   YS <- YS[!badRow]
   XO <- XO[!badRow,]
   YO <- YO[!badRow]
   probitDummy <- probitDummy[!badRow]
   ## Now indices for packing the separate outcomes into full outcome vectors.  Note we treat
   ## invMillsRatio as a separate parameter
   iBetaS <- seq(length=NXS)
   iBetaO <- NXS + seq(length=NXO)
   iMills <- NXS + NXO + 1
   iSigma <- iMills + 1
   iRho <- iSigma + 1
   ##
   nObs <- length(YS)
   nParam <- iRho
   N0 <- sum(YS == levels(YS)[1])
   N1 <- nObs - N0
                                        # sigma, rho
   if( print.level > 0 ) {
      cat ( "\nEstimating 1st step Probit model . . ." )
   }
   result$probit <- probit(YS ~ XS - 1, x=TRUE, print.level=print.level - 1, iterlim=30)
                                        # a large iterlim may help with weakly identified models
   if( print.level > 0 ) {
       cat( " OK\n" )
   }
   imrData <- invMillsRatio( result$probit )
   result$imrDelta <- drop(imrData$delta1)
   ## ---- Outcome estimation -----
   if( is.null( inst ) ) {
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step (outcome) OLS model . . ." )
      }
      outcomeMod <- lm(YO ~ -1 + XO + imrData$IMR1,
                      subset = probitDummy )
      intercept <- any(apply(model.matrix(outcomeMod), 2,
                             function(v) (v[1] > 0) & (all(v == v[1]))))
                                        # we have determine whether the outcome model has intercept.
                                        # This is necessary later for calculating R^2
      resid <- residuals( outcomeMod )
      step2coef <- coef( outcomeMod )
      names(step2coef) <- c(colnames(XO), "invMillsRatio")
      if(print.level > 0)
          cat( " OK\n" )
   }
   else {
      data$invMillsRatio <- imrData$IMR1
      if( print.level > 0 ) {
         cat ( "Estimating 2nd step 2SLS/IV model . . ." )
      }
      step2formula <- as.formula( paste( outcome[ 2 ], "~", outcome[ 3 ],
                                        "+ invMillsRatio" ) )
      formulaList <- list( step2formula )
      instImr <- as.formula( paste( "~", inst[ 2 ], "+ invMillsRatio" ) )
      library( systemfit )
      outcomeMod <- systemfit( "2SLS", formulaList, inst = instImr,
                             data = data[ probitDummy, ] )
      intercept = FALSE
                                        # we calculate R^2 differently here (hopefully)
      resid <- residuals( outcomeMod )[ , 1 ]
      step2coef <- coefficients( outcomeMod$eq[[ 1 ]] )
      if( print.level > 0 ) cat( " OK\n" )
   }
   result$sigma <- drop( sqrt( crossprod( resid ) /
                              sum( probitDummy ) +
                              mean(result$imrDelta[ probitDummy ] ) *
                              step2coef[ "invMillsRatio" ]^2 ) )
   result$rho <-  step2coef[ "invMillsRatio" ] / result$sigma
   names(result$rho) <- NULL
   result$invMillsRatio <- drop( imrData$IMR1 )
   ## Stack all final coefficients to 'coefficients'
   coefficients <- c(coef(result$probit),
                     step2coef[names(step2coef) != "invMillsRatio"],
                     step2coef["invMillsRatio"],
                     sigma=result$sigma, rho=result$rho)
   names(coefficients)[iBetaS] <- gsub("^XS", "", names(coefficients)[iBetaS])
   if( print.level > 0 ) {
      cat ( "Calculating coefficient covariance matrix . . ." )
   }
                                        # the following variables are named according to Greene (2003), p. 785
   if( is.null( inst ) ) {
      xMat <- model.matrix( outcomeMod )
   } else {
      xMat <- outcomeMod$eq[[ 1 ]]$x
   }
   ## Varcovar matrix.  Fill only a few parts, rest will remain NA
   vc <- matrix(0, nParam, nParam)
   colnames(vc) <- row.names(vc) <- names(coefficients)
   vc[] <- NA
   if(!is.null(vcov(result$probit)))
       vc[iBetaS,iBetaS] <- vcov(result$probit)
   vc[c(iBetaO, iMills), c(iBetaO, iMills)] <- heckitVcov( xMat,
                                                          model.matrix( result$probit )[ probitDummy, ],
                                                          vcov( result$probit ),
                                                          result$rho,
                                                          result$imrDelta[ probitDummy ],
                                                          result$sigma )
                                        # here we drop invMillsRatio part
   result$vcov <- vc
   ##
   if( print.level > 0 )
       cat( " OK\n" )
   result$coefficients <- coefficients
                                        # for coef() etc. methods
   ## the 'param' component is intended to all kind of technical info
   result$param <- list(index=list(betaS=iBetaS, betaO=iBetaO, 
                        Mills=iMills, sigma=iSigma, rho=iRho,
                        errTerms = c( iMills, iSigma, iRho ),
                        outcome = c( iBetaO, iMills ) ),
                                        # The location of results in the coef vector
                        oIntercept=intercept,
                        N0=N0, N1=N1,
                        nParam=nParam, nObs=nObs, df=nObs-nParam+1)
   result$lm <- outcomeMod
   result$tobitType <- 2
   result$method <- "2step"
   class( result ) <- c( "selection", class(result))
   return( result )
}
