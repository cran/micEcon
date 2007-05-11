tobit2 <- function(selection, formula,
                   data=sys.frame(sys.parent()),
                   method="ml",
                   start=NULL, print.level=0,
                   y1=FALSE, z=FALSE, y2=FALSE, x=FALSE, model=FALSE,
                   ...) {
   ## The model (Amemiya 1985):
   ## The latent variables are:
   ## y1* = z'gamma + u1
   ## y2* = x'beta + u2
   ## The observables are:
   ##      / 1  if  y1* > 0
   ## y1 = \ 0  if  y1* <= 0
   ##      / 0    if  y1 = 0
   ## y2 = \ y2*  if  y1 = 1
   ##
   ## PARAMETERS:
   ## selection: formula for the selection equation (y1)
   ## formula    main model
   ## method     maximum likelihood (ml) or two-step (2step)
   ## data       dataset
   ## start      initial value of coefficients.  The order is as follows:
   ##            start = (gamma', beta', sigma, rho)'
   ## print.level: 0 - nothing printed, as larger, as more information printed
   ## y1..z        whether to return corresponding matrixis of explanatory and dependent variables
   ## model        whether to return model frames
   ##  ...         additional parameters form maximisation
   ##
   ## RESULTS:
   ## a list with the following components:
   ##
   loglik <- function( beta) {
      g <- beta[iGamma]
      b <- beta[iBeta]
      sigma <- beta[iSigma]
      if(sigma < 0) return(NA)
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
      Z1.g <- Z1 %*% g
      Z2.g <- Z2 %*% g
      X2.b <- X2 %*% b
      u2 <- Y21 - X2.b
      r <- sqrt( 1 - rho^2)
      B <- (Z2.g + rho/sigma*u2)/r
      l1 <- sum(pnorm(-Z1.g, log.p=TRUE))
                                        # loglik related with unobserved cases ...
      l2 <- -N2/2*log(2*pi) - N2*log(sigma) +
         sum(pnorm(B, log.p=TRUE) - 0.5*(u2/sigma)^2)
                                        # ... and observed cases
      loglik <- l1 + l2
   }
   gradlik <- function(beta) {
      g <- beta[iGamma]
      b <- beta[iBeta]
      sigma <- beta[iSigma]
      if(sigma < 0) return(NA)
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
      Z1.g <- Z1 %*% g
      Z2.g <- Z2 %*% g
      X2.b <- X2 %*% b
      u2 <- Y21 - X2.b
      r <- sqrt( 1 - rho^2)
      B <- (Z2.g + rho/sigma*u2)/r
      lambdaB <- dnorm(B)/pnorm(B)
      gradient <- numeric(nParam)
      gradient[iGamma] <- t(Z1) %*% (-dnorm(-Z1.g)/pnorm(-Z1.g)) +
         (t(Z2) %*% lambdaB)/r
      gradient[iBeta] <- t(X2) %*% (u2/sigma^2 - lambdaB*rho/sigma/r)
      gradient[iSigma] <- sum(u2^2/sigma^3 - lambdaB*rho*u2/sigma^2/r) - N2/sigma
      gradient[iRho] <- sum(lambdaB*(u2/sigma + rho*Z2.g))/r^3
      gradient
   }
   hesslik <- function(beta) {
      g <- beta[iGamma]
      b <- beta[iBeta]
      sigma <- beta[iSigma]
      if(sigma < 0) return(NA)
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
      Z1.g <- as.vector(Z1 %*% g)
      Z2.g <- as.vector(Z2 %*% g)
      X2.b <- as.vector(X2 %*% b)
      u2 <- Y21 - X2.b
      r <- sqrt( 1 - rho^2)
      B <- (Z2.g + rho/sigma*u2)/r
      lambdaB <- dnorm(B)/pnorm(B)
      fZ1.g <- dnorm(-Z1.g)
      FZ1.g <- pnorm(-Z1.g)
                                        # the previous code
                                        #       C <- ifelse(B > -25,
                                        #                   -(pnorm(B)*dnorm(B)*B + dnorm(B)^2)/pnorm(B)^2,
                                        #                   -1)
      C <- ifelse(B > -500,
                  -exp(dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))*B -
                  exp(2 * (dnorm(B, log = TRUE) - pnorm(B, log.p = TRUE))),
                  -1)
                                        # recommended by Dimitrios Rizopoulos, KULeuven
                                        # This is a hack in order to avoid numerical problems.  How to do
                                        # it better?  How to prove the limit value?
      hess <- matrix(0, nParam, nParam)
      a <- (-fZ1.g*FZ1.g*Z1.g + fZ1.g^2)/FZ1.g^2
      hess[iGamma,iGamma] <- -t(Z1) %*% (Z1*a) + t(Z2) %*% (Z2*C)/r^2
      hess[iGamma,iBeta] <- -t(Z2) %*% (X2*C)*rho/r^2/sigma
      hess[iBeta,iGamma] <- t(hess[iGamma,iBeta])
      hess[iGamma,iSigma] <- -rho/sigma^2/r^2*t(Z2) %*% (C*u2)
      hess[iSigma,iGamma] <- t(hess[iGamma,iSigma])
      hess[iGamma,iRho] <- t(Z2) %*%
         (C*(u2/sigma + rho*Z2.g)/r^4 + lambdaB*rho/r^3)
      hess[iRho,iGamma] <- t(hess[iGamma,iRho])
      hess[iBeta,iBeta] <- t(X2) %*%
         (X2 * ((rho/r)^2*C - 1))/sigma^2
      hess[iBeta,iSigma] <- t(X2) %*%
         (C*rho^2/sigma^3*u2/r^2 +
            rho/sigma^2*lambdaB/r - 2*u2/sigma^3)
      hess[iSigma,iBeta] <- t(hess[iBeta,iSigma])
      hess[iBeta,iRho] <- t(X2) %*%
         (-C*(u2/sigma + rho*Z2.g)/r^4*rho -
            lambdaB/r^3)/sigma
      hess[iRho,iBeta] <- t(hess[iBeta,iRho])
      hess[iSigma,iSigma] <- sum(
                                   -3*u2*u2/sigma^4
                                   +2*lambdaB* u2/r *rho/sigma^3
                                   +rho^2/sigma^4 *u2*u2/r^2 *C) +
                                       N2/sigma^2
      hess[iSigma,iRho] <- hess[iRho,iSigma] <-
         -sum((C*rho*(u2/sigma + rho*Z2.g)/r + lambdaB)*
         u2/sigma^2)/r^3
      hess[iRho,iRho] <-
         sum(C*((u2/sigma + rho*Z2.g)/r^3)^2 +
         lambdaB*(Z2.g*(1 + 2*rho^2) + 3*rho*u2/sigma) / r^5 )
      return( hess )
   }
   ## --- the main program ---
   ## First the consistency checks
   .Deprecated( "selection", "micEcon" )
   if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   }
   if( length( formula ) != 3 ) {
      stop( "argument 'formula' must be a 2-sided formula" )
   }
   if( "probit" %in% substr( all.vars( formula ), 1, 6 ) ) {
      stop( "argument 'formula' may not include variable names",
                  " starting with 'probit'" )
   }
   if( class( selection ) != "formula" ) {
      stop( "argument 'selection' must be a formula" )
   }
   if( length( selection ) != 3 ) {
      stop( "argument 'selection' must be a 2-sided formula" )
   }
   if( "probit" %in% substr( all.vars( selection ), 1, 6 ) ) {
      stop( "argument 'selection' may not include a variable",
                  " names starting with 'probit'" )
   }
   probitEndogenous <- model.frame( selection, data = data )[ , 1 ]
   probitLevels <- levels( as.factor( probitEndogenous ) )
   if( length( probitLevels ) != 2 ) {
      stop( "the left hand side of 'selection' has to contain",
         " exactly two levels (e.g. FALSE and TRUE)" )
   }
   data$probitDummy <- probitEndogenous == probitLevels[ 2 ]
   ## now check whether two-step method is needed: either for final estimate or initial parameters
   if(method == "2step" | is.null(start)) {
      twoStep <- heckit2fit(selection, formula, data)
      if(method == "2step") {
         return(twoStep)
      }
   }
   ## Now extract model frames etc
   ## Y1 (selection equation)
   cl <- match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("selection", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mf1 <- mf[c(1, m)]
   mf1$drop.unused.levels <- TRUE
   mf1[[1]] <- as.name("model.frame")
   names(mf1)[2] <- "formula"
                                        # model.frame requires the parameter to be 'formula'
   mf1 <- eval(mf1, parent.frame())
   mt1 <- attr(mf1, "terms")
   Z <- model.matrix(mt1, mf1)
   Y1 <- model.response(mf1, "numeric")
   ## Y2 (regression)
   m <- match(c("formula", "data", "subset", "weights", "na.action",
                "offset"), names(mf), 0)
   mf2 <- mf[c(1, m)]
   mf2$drop.unused.levels <- TRUE
   mf2[[1]] <- as.name("model.frame")
   mf2 <- eval(mf2, parent.frame())
                                        # Note: if unobserved variables are marked as NA, eval returns a
                                        # subframe of visible variables only.  We have to check it later
   mt2 <- attr(mf2, "terms")
   X <- model.matrix(mt2, mf2)
   Y2 <- model.response(mf2, "numeric")
                                        #
   NZ <- ncol( Z)
   NX <- ncol( X)
   nParam <- NZ + NX + 2
                                        # Total # of parameters
   nObs <- length(Y1)
   N1 <- length( Y1[Y1==0])
   N2 <- length( Y1[Y1==1])
   ## indices in for the parameter vector
   iGamma <- 1:NZ
   iBeta <- max(iGamma) + seq(length=NX)
   iSigma <- max(iBeta) + 1
   iRho <- max(iSigma) + 1
   ## divide data by choices
   if(dim(X)[1] == N2) {
                                        # unobserved variables marked as NA
      Y21 <- Y2
      X2 <- X
    }
   else {
                                        # unobserved variables are (technically) correct numbers, e.g.
                                        # zeros
      Y21 <- Y2[Y1==1]
      X2 <- X[Y1==1,,drop=FALSE]
   }
   Z1 <- Z[Y1==0,,drop=FALSE]
   Z2 <- Z[Y1==1,,drop=FALSE]
   if(print.level > 0) {
      cat( "Not observed:", N1, "; observed:", N2,"\n", sep="")
   }
   ## initial values for parameters.  Note that 'twoStep' is calculated before
   if(is.null(start)) {
      if(print.level > 0) {
         cat("Initial values by 2-step method:results\n")
         print(twoStep)
      }
      start[iGamma] <- coef(twoStep)[twoStep$param$index$betaS]
      start[iBeta] <- coef(twoStep)[twoStep$param$index$betaO]
      start[iSigma] <- coef(twoStep)[twoStep$param$index$sigma]
      start[iRho] <- coef(twoStep)[twoStep$param$index$rho]
      if(start[iRho] > 0.99)
        start[iRho] <- 0.99
      else if(start[iRho] < -0.99)
        start[iRho] <- -0.99
      names(start) <- c(colnames(Z), colnames(X), "sigma", "rho")
      if(print.level > 0) {
         cat("Initial values:\n")
         print(start)
      }
   }
   estimation <- maxLik(loglik, gradlik, hesslik,
                        start=start,
                        print.level=print.level - 1, ...)
   ## The following commented lines are for testing analytic gradient and Hessian
#    compare.derivatives(loglik, gradlik, t0=init, ...)
#    compare.derivatives(gradlik, hesslik, t0=init, ...)
   param <- list(index=list(betaS=iGamma, betaO=iBeta,
                   sigma=iSigma, rho=iRho),
                 nParam=nParam,
                 nObs=nObs,
                 N0=N1,
                 N1=N2,
                 NXS=NZ,
                 NXO=NX,
                 df=nObs - estimation$NActiveParam)
   result <- c(estimation,
               twoStep=list(twoStep),
               param=list(param),
               call=cl,
               terms1=mt1,
               terms2=mt2,
               y1=switch(y1, "1"=list(Y1), "0"=NULL),
               z=switch(z, "1"=list(X), "0"=NULL),
               y2=switch(y2, "1"=list(Y2), "0"=NULL),
               x=switch(x, "1"=list(X), "0"=NULL),
               model=switch(model, "1"=list(selection=mf1, formula=mf2), "0"=NULL))
   result$tobitType <- 2
   class(result) <- c("selection", class(estimation))
   return(result)
}
