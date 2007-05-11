probit <- function( formula, subset, start=NULL,
                   data=sys.frame(sys.parent()),
                   x=FALSE, y=FALSE, model=FALSE,
                   method="ML", 
                   ...) {
   ## Probit binary choice model
   ## formula: model formula, response must be either a logical or numeric vector containing only 0-s and
   ##          1-s
   ## start:      initial value of the parameters
   ## data     dataframe where the variables are defined
   ## x        whether to return model matrix
   ## y                                response vector
   ## model                            frame
   ## method   method for evaluation:
   ##          ML            maximum likelihood
   ##          model.frame   don't evaluate, only return the frame
   ## ...      further arguments for the maxLik algorithm
   ##
   ## return: a list with following components.
   ##  $results: maximisation results
   ##  $LRT:     list with two components:
   ##            LRT: likelihood ration test H0 - none of the variables significant
   ##            df:  corresponding degrees of freedom
   ##  x         model matrix (only if requested)
   ##  call      call
   ##  terms     terms
   loglik <- function( beta) {
      xb0 <- x0 %*% beta
      xb1 <- x1 %*% beta
      loglik <- sum(pnorm( xb0, log=TRUE, lower.tail=FALSE)) + sum(pnorm( xb1, log.p=TRUE))
   }
   gradlik <- function(beta) {
      ## gradient is 1 x nParam matrix
      xb0 <- x0 %*% beta
      xb1 <- x1 %*% beta
      gradlik <- - t(dnorm(xb0)/pnorm(xb0, lower.tail=FALSE)) %*% x0 +
          t(dnorm(xb1)/pnorm(xb1)) %*% x1
   }
   hesslik <- function(beta) {
      xb0 <- as.vector( x0 %*% beta)
      xb1 <- as.vector( x1 %*% beta)
      F0 <- pnorm( xb0)
      F1 <- pnorm( xb1)
      f0 <- dnorm( xb0)
      f1 <- dnorm( xb1)
      yF0 <- pnorm( xb0, lower.tail=FALSE)
      loglik2 <- t( x1) %*% ( x1 * ( -f1*xb1*F1 - f1*f1)/F1/F1) -
         t( x0) %*% ( x0 * ( -f0*xb0*yF0 + f0*f0)/yF0/yF0)
                    # note that df/db' = -f (x'b) x'
   }
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf$drop.unused.levels <- TRUE
   mf[[1]] <- as.name("model.frame")
   eval(data)
                                        # we need to eval() data here, otherwise the evaluation of the
                                        # model frame will be wrong if called from inside a function
                                        # inside a function (sorry, I can't understand it either :-(
   mf <- eval(mf, envir=parent.frame())
   if (method == "model.frame")
       return(mf)
   else if (method != "ML")
       warning("method = ", method, " is not supported. Using \"ML\"")
   mt <- attr(mf, "terms")
   Y <- model.response( mf )
   X <- model.matrix(mt, mf, contrasts)
   nParam <- ncol( X)
   nObs <- length( Y)
   N1 <- length( y[Y != 0])
   N0 <- nObs - N1
   if(N0 == 0 | N1 == 0) {
      stop("No variance in the response variable")
   }
   x0 <- X[Y==0,]
   x1 <- X[Y==1,]
   if(is.null(start)) {
      start <- rep( 0, nParam)
   }
   if(is.null(names(start))) {
      names(start) <- dimnames(X)[[2]]
   }
   ## Main estimation
   estimation <- maxLik(loglik, gradlik, hesslik, start,
                        method="Newton-Raphson", ...)
   ## compare.derivatives(gradlik, hesslik, t0=start)
                                        #
   ## Likelihood ratio test: H0 -- all the coefficients, except intercept
   ## are zeros.  ML estimate for this model is qnorm(N1/nObs)
   ll.bar <- loglik(c(qnorm(N1/nObs), rep(0, nParam-1)))
   LRT <- 2*(estimation$maximum - ll.bar)
                                        #
   result <- c(estimation,
               LRT=list(list(LRT=LRT, df=nParam-1)),
                                        # there are df-1 constraints
               param=list(list(nParam=nParam,nObs=nObs, N1=N1, N0=N0)),
               df=nObs - nParam,
               call=cl,
               terms=mt,
               x=switch(x, "1"=list(X), "0"=NULL),
               y=switch(y, "1"=list(Y), "0"=NULL),
               model=switch(model, "1"=list(mf), "0"=NULL))
   class(result) <- c("probit", class(estimation))
   result
}
