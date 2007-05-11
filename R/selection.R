selection <- function(selection, outcome,
                      data=sys.frame(sys.parent()),
                      method="ml",
                      start=NULL,
                      ys=FALSE, xs=FALSE,
                      yo=FALSE, xo=FALSE,
                      mfs=FALSE, mfo=FALSE,
                      print.level=0,
                      ...) {
   ## Heckman-style sample-selection models
   ## selection:   formula
   ##              LHS: must be convertable to two-level factor (e.g. 0-1, 1-2, "A"-"B")
   ##              RHS: ordinary formula as in lm()
   ## outcome:     formula, or a list of two formulas.
   ##              If the outcome contains one formula, it is
   ##              assumed that we observe the formula only if selection equals to the first level
   ##              of selection (e.g. "0", "1" or "A" in the examples above).
   ##              If the outcome contains two formulas, we assume that the first one is observed if
   ##              selection equals to the first level, otherwise the second formula.
   ## ys, xs, yo, xo, mfs, mfo: whether to return the response, model matrix or
   ##              the model frame of outcome and selection equation(s)
   ## First the consistency checks
   ## ...          additional arguments for tobit2fit and tobit5fit
   type <- 0
   if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
                                        #
   if("formula" %in% class( outcome )) {
      if( length( outcome ) != 3 ) {
         stop( "argument 'outcome' must be a 2-sided formula" )
      }
      type <- 2
   }
   else if("list" %in% class(outcome)) {
      if(length(outcome) == 1) {
         outcome <- outcome[[1]]
         type <- 2
      }
      else if(length(outcome) == 2) {
         if("formula" %in% class(outcome[[1]])) {
            if( length( outcome[[1]] ) != 3 ) {
               stop( "argument 'outcome[[1]]' must be a 2-sided formula" )
            }
         }
         else
             stop( "argument 'outcome[[1]]' must be either a formula or a list of two formulas" )
         if("formula" %in% class(outcome[[2]])) {
            if( length( outcome[[2]] ) != 3 ) {
               stop( "argument 'outcome[[2]]' must be a 2-sided formula" )
            }
         }
         else
             stop( "argument 'outcome[[2]]' must be either a formula or a list of two formulas" )
         type <- 5
      }
      else
          stop("argument 'outcome' must contain 1 or 2 components")
   }
   else
       stop( "argument 'selection' must be either a formula or a list of two formulas" )
   if(print.level > 0)
       cat("Tobit", type, "model\n")
   probitEndogenous <- model.frame( selection, data = data )[ , 1 ]
   probitLevels <- levels( as.factor( probitEndogenous ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   data$probitDummy <- probitEndogenous == probitLevels[ 2 ]
   ## now check whether two-step method was requested
   cl <- match.call()
   if(method == "2step") {
      if(type == 2)
          twoStep <- heckit2fit(selection, outcome, data=data,
            print.level = print.level, ... )
      else if(type == 5)
          twoStep <- heckit5fit(selection, outcome, data=data,
            print.level = print.level, ... )
      else
          stop("unknown type")
      twoStep$call <- cl
      class(twoStep) <- c("selection", class(twoStep))
      return(twoStep)
   }
   ## Now extract model frames etc
   ## YS (selection equation)
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset"), names(mf), 0)
   mfS <- mf[c(1, m)]
   mfS$drop.unused.levels <- TRUE
   mfS$na.action <- na.pass
   mfS[[1]] <- as.name("model.frame")
   names(mfS)[2] <- "formula"
                                        # model.frame requires the parameter to
                                        # be 'formula'
   mfS <- eval(mfS, parent.frame())
   mtS <- attr(mfS, "terms")
   XS <- model.matrix(mtS, mfS)
   YS <- model.response(mfS, "numeric")
   ## check for NA-s.  Because we have to find NA-s in several frames, we cannot use the standard na.
   ## functions here.  Find bad rows and remove them later.
   ## We check XS and YS separately, because mfS may be a data frame with complex structure (e.g.
   ## including matrices)
   badRow <- is.na(YS)
   badRow <- badRow | apply(XS, 1, function(v) any(is.na(v)))
   ## YO (outcome equation)
   if(type == 2) {
      oArg <- match("outcome", names(mf), 0)
                                        # find the outcome argument
      m <- match(c("outcome", "data", "subset", "weights", 
                   "offset"), names(mf), 0)
      ## replace the outcome list by the first equation and evaluate it
      mfO <- mf[c(1, m)]
      mfO$drop.unused.levels <- TRUE
      mfO$na.action <- na.pass
      mfO[[1]] <- as.name("model.frame")
                                        # eval it as model frame
      names(mfO)[2] <- "formula"
      mfO <- eval(mfO, parent.frame())
                                        # Note: if unobserved variables are
                                        # marked as NA, eval returns a
                                        # subframe of visible variables only.
                                        # We have to check it later
      mtO <- attr(mfO, "terms")
      XO <- model.matrix(mtO, mfO)
      YO <- model.response(mfO, "numeric")
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
      NXS <- ncol(XS)
      NXO <- ncol(XO)
      iGamma <- 1:NXS
      iBeta <- max(iGamma) + seq(length=NXO)
      iSigma <- max(iBeta) + 1
      iRho <- max(iSigma) + 1
      twoStep <- NULL
      if(is.null(start)) {
         twoStep <- heckit2fit(selection, outcome, data=data,
            print.level = print.level, ... )
         coefs <- coef(twoStep, part="full")
         start[iGamma] <- coefs[twoStep$param$index$betaS]
         start[iBeta] <- coefs[twoStep$param$index$betaO]
         start[iSigma] <- coefs[twoStep$param$index$sigma]
         start[iRho] <- coefs[twoStep$param$index$rho]
         if(start[iRho] > 0.99)
             start[iRho] <- 0.99
         else if(start[iRho] < -0.99)
             start[iRho] <- -0.99
         names(start) <- c(colnames(XS), colnames(XO), "sigma", "rho")
      }
      estimation <- tobit2fit(YS, XS, YO, XO, start,
                              print.level=print.level, ...)
      param <- list(index=list(betaS=iGamma,
                    betaO=iBeta, sigma=iSigma, rho=iRho,
                    errTerms = c( iSigma, iRho ),
                    outcome = iBeta ),
                    NXS=ncol(XS), NXO=ncol(XO),
                    N0=sum(YS==0), N1=sum(YS==1),
                    nObs=length(YS), nParam=length(start),
                    df=length(YS) - length(start))
   }
   else if(type == 5) {
      oArg <- match("outcome", names(mf), 0)
                                        # find the outcome argument
      ocome <- as.list(mf[[oArg]])
      formula1 <- ocome[[2]]
      formula2 <- ocome[[3]]
                                        # Now we have extracted both formulas
      m <- match(c("outcome", "data", "subset", "weights",
                   "offset"), names(mf), 0)
      ## replace the outcome list by the first equation and evaluate it
      mf[[oArg]] <- formula1
      mf1 <- mf[c(1, m)]
      mf1$drop.unused.levels <- TRUE
      mf1$na.action = na.pass
      mf1[[1]] <- as.name("model.frame")
                                        # eval it as model frame
      names(mf1)[2] <- "formula"
      mf1 <- eval(mf1, parent.frame())
      mtO1 <- attr(mf1, "terms")
      XO1 <- model.matrix(mtO1, mf1)
      YO1 <- model.response(mf1, "numeric")
      badRow <- badRow | (is.na(YO1) & (!is.na(YS) & YS == 0))
      badRow <- badRow | (apply(XO1, 1, function(v) any(is.na(v))) & (!is.na(YS) & YS == 0))
      ## repeat all the stuff with second equation
      mf[[oArg]] <- formula2
      mf2 <- mf[c(1, m)]
      mf2$drop.unused.levels <- TRUE
      mf2$na.action <- na.pass
      mf2[[1]] <- as.name("model.frame")
                                        # eval it as model frame
      names(mf2)[2] <- "formula"
      mf2 <- eval(mf2, parent.frame())
      mtO2 <- attr(mf2, "terms")
      XO2 <- model.matrix(mtO2, mf2)
      YO2 <- model.response(mf2, "numeric")
      badRow <- badRow | (is.na(YO2) & (!is.na(YS) & YS == 1))
      badRow <- badRow | (apply(XO2, 1, function(v) any(is.na(v))) & (!is.na(YS) & YS == 1))
      ## indices in for the parameter vector.  These are returned in order to provide the user a way
      ## to extract certain components from the coefficients
      NXS <- ncol(XS)
      NXO1 <- ncol(XO1)
      NXO2 <- ncol(XO2)
      XS <- XS[!badRow,]
      YS <- YS[!badRow]
      XO1 <- XO1[!badRow,]
      YO1 <- YO1[!badRow]
      XO2 <- XO2[!badRow,]
      YO2 <- YO2[!badRow]
      iBetaS <- 1:NXS
      iBetaO1 <- seq(tail(iBetaS, 1)+1, length=NXO1)
      iSigma1 <- tail(iBetaO1, 1) + 1
      iRho1 <- tail(iSigma1, 1) + 1
      iBetaO2 <- seq(tail(iRho1, 1) + 1, length=NXO2)
      iSigma2 <- tail(iBetaO2, 1) + 1
      iRho2 <- tail(iSigma2, 1) + 1
      nParam <- iRho2
      twoStep <- NULL
      if(is.null(start)) {
         start <- numeric(nParam)
         if(print.level > 0) {
            cat("Start values by Heckman 2-step method (", nParam, " componenets)\n", sep="")
         }
         twoStep <- heckit5fit(selection, as.formula(formula1), as.formula(formula2),
                           data=data, print.level = print.level, ... )
         ind <- twoStep$param$index
         start <- coef(twoStep, part="full")[c(ind$betaS,
                                  ind$betaO1, ind$sigma1, ind$rho1,
                                  ind$betaO2, ind$sigma2, ind$rho2)]
         names( start ) <- sub( "^[SO][12]?:", "", names( start ) )
      }
      estimation <- tobit5fit(YS, XS, YO1, XO1, YO2, XO2, start=start,
                              print.level=print.level, ...)
      param <- list(index=list(betaS=iBetaS,
                    betaO1=iBetaO1, sigma1=iSigma1, rho1=iRho1,
                    betaO2=iBetaO2, sigma2=iSigma2, rho2=iRho2,
                    errTerms = c( iSigma1, iSigma2, iRho1, iRho2 ),
                    outcome = c( iBetaO1, iBetaO2 ) ),
                    NXS=ncol(XS),
                    NXO1=ncol(XO1), NXO2=ncol(XO2),
                    N1=sum(YS==0), N2=sum(YS==1),
                    nObs=length(YS), nParam=length(start),
                    df=length(YS) - length(start))
   }
   ## now fit the model
   result <- c(estimation,
               twoStep=list(twoStep),
               start=list(start),
               param=list(param),
               call=cl,
               termsS=mtS,
               termsO=switch(as.character(type), "1"=mtO, "5"=list(mtO1, mtO2), "0"=NULL),
               ys=switch(as.character(ys), "TRUE"=list(YS), "FALSE"=NULL),
               xs=switch(as.character(xs), "TRUE"=list(XS), "FALSE"=NULL),
               yo=switch(as.character(yo),
               "TRUE"=switch(as.character(type), "2"=list(YO), "5"=list(YO1,
                                                               YO2)), "FALSE"=NULL),
               xo=switch(as.character(xo),
               "TRUE"=switch(as.character(type), "2"=list(XO), "5"=list(XO1, XO2)), "FALSE"=NULL),
               mfs=switch(as.character(mfs), "TRUE"=list(mfS), "FALSE"=NULL),
               mfo=switch(as.character(mfs),
               "TRUE"=switch(as.character(type), "2"=list(mfO), "5"=list(mf1,
                      mf2), "FALSE"=NULL))
               )
   class( result ) <- class( estimation ) 
   return(result)
}
