tobit2fit <- function(YS, XS, YO, XO, start,
                      print.level=0,
                      maxMethod="Newton-Raphson",
                      ...) {
### The model is as follows (Amemiya 1985):
### The latent variables are:
### YS* = XS'g + u1
### y2* = X'b2 + u2
### y3* = X'b3 + u3
### The observables are:
###      / 1  if  YS* > 0
### YS = \ 0  if  YS* <= 0
###      / y2*  if  YS = 0
### y2 = \ 0    if  YS = 1
###      / 0    if  YS = 0
### y3 = \ y3*  if  YS = 1
### 
###  YS      binary or logical vector, 0 (FALSE) corresponds to observable y2*,
###          1 (TRUE) to observable y3*
###  y2, y3  numeric vector, outcomes
###  X       explanatory variables for outcomes
###  XS              -"-                selection, should include exclusion restriction
###  ...     additional parameters for maxLik
### Result:
### Object of class 'tobit2', derived from 'maxLik'.
### Includes all the components of maxLik and additionally
### twoStep   Results for Heckman two-step estimator, used for initial values
### 
   loglik <- function( beta) {
      g <- beta[iBetaS]
      b <- beta[iBetaO]
      sigma <- beta[iSigma]
      if(sigma < 0) return(NA)
      rho <- beta[iRho]
      if( ( rho < -1) || ( rho > 1)) return(NA)
      XS0.g <- XS0 %*% g
      XS1.g <- XS1 %*% g
      XO1.b <- XO1 %*% b
      u2 <- YO1 - XO1.b
      r <- sqrt( 1 - rho^2)
      B <- (XS1.g + rho/sigma*u2)/r
      l1 <- sum(pnorm(-XS0.g, log.p=TRUE))
                                        # loglik related with unobserved cases ...
      l2 <- -N1/2*log(2*pi) - N1*log(sigma) +
          sum(pnorm(B, log.p=TRUE) - 0.5*(u2/sigma)^2)
                                        # ... and observed cases
      loglik <- l1 + l2
   }
    gradlik <- function(beta) {
       g <- beta[iBetaS]
       b <- beta[iBetaO]
       sigma <- beta[iSigma]
       if(sigma < 0) return(NA)
       rho <- beta[iRho]
       if( ( rho < -1) || ( rho > 1)) return(NA)
       XS0.g <- XS0 %*% g
       XS1.g <- XS1 %*% g
       XO1.b <- XO1 %*% b
       u2 <- YO1 - XO1.b
       r <- sqrt( 1 - rho^2)
       B <- (XS1.g + rho/sigma*u2)/r
       lambdaB <- ifelse(B > -30, dnorm(B)/pnorm(B), B)
                                        # This is a hack in order to avoid numeric problems
       gradient <- numeric(nParam)
       gradient[iBetaS] <- t(XS0) %*% (-dnorm(-XS0.g)/pnorm(-XS0.g)) +
           (t(XS1) %*% lambdaB)/r
       gradient[iBetaO] <- t(XO1) %*% (u2/sigma^2 - lambdaB*rho/sigma/r)
       gradient[iSigma] <- sum(u2^2/sigma^3 - lambdaB*rho*u2/sigma^2/r) - N1/sigma
       gradient[iRho] <- sum(lambdaB*(u2/sigma + rho*XS1.g))/r^3
       gradient
    }
    hesslik <- function(beta) {
       g <- beta[iBetaS]
       b <- beta[iBetaO]
       sigma <- beta[iSigma]
       if(sigma < 0) return(NA)
       rho <- beta[iRho]
       if( ( rho < -1) || ( rho > 1)) return(NA)
       XS0.g <- as.vector(XS0 %*% g)
       XS1.g <- as.vector(XS1 %*% g)
       XO1.b <- as.vector(XO1 %*% b)
       u2 <- YO1 - XO1.b
       r <- sqrt( 1 - rho^2)
       B <- (XS1.g + rho/sigma*u2)/r
       lambdaB <- ifelse(B > -30, dnorm(B)/pnorm(B), B)
                                        # This is a hack in order to avoid numeric problems
       fXS0.g <- dnorm(-XS0.g)
       FXS0.g <- pnorm(-XS0.g)
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
       a <- (-fXS0.g*FXS0.g*XS0.g + fXS0.g^2)/FXS0.g^2
       hess[iBetaS,iBetaS] <- -t(XS0) %*% (XS0*a) + t(XS1) %*% (XS1*C)/r^2
       hess[iBetaS,iBetaO] <- -t(XS1) %*% (XO1*C)*rho/r^2/sigma
       hess[iBetaO,iBetaS] <- t(hess[iBetaS,iBetaO])
       hess[iBetaS,iSigma] <- -rho/sigma^2/r^2*t(XS1) %*% (C*u2)
       hess[iSigma,iBetaS] <- t(hess[iBetaS,iSigma])
       hess[iBetaS,iRho] <- t(XS1) %*%
           (C*(u2/sigma + rho*XS1.g)/r^4 + lambdaB*rho/r^3)
       hess[iRho,iBetaS] <- t(hess[iBetaS,iRho])
       hess[iBetaO,iBetaO] <- t(XO1) %*%
           (XO1 * ((rho/r)^2*C - 1))/sigma^2
       hess[iBetaO,iSigma] <- t(XO1) %*%
           (C*rho^2/sigma^3*u2/r^2 +
            rho/sigma^2*lambdaB/r - 2*u2/sigma^3)
       hess[iSigma,iBetaO] <- t(hess[iBetaO,iSigma])
       hess[iBetaO,iRho] <- t(XO1) %*%
           (-C*(u2/sigma + rho*XS1.g)/r^4*rho -
            lambdaB/r^3)/sigma
       hess[iRho,iBetaO] <- t(hess[iBetaO,iRho])
       hess[iSigma,iSigma] <- sum(
                                  -3*u2*u2/sigma^4
                                  +2*lambdaB* u2/r *rho/sigma^3
                                  +rho^2/sigma^4 *u2*u2/r^2 *C) +
                                      N1/sigma^2
       hess[iSigma,iRho] <- hess[iRho,iSigma] <-
           -sum((C*rho*(u2/sigma + rho*XS1.g)/r + lambdaB)*
                u2/sigma^2)/r^3
       hess[iRho,iRho] <-
           sum(C*((u2/sigma + rho*XS1.g)/r^3)^2 +
               lambdaB*(XS1.g*(1 + 2*rho^2) + 3*rho*u2/sigma) / r^5 )
       return( hess )
    }
    ## ---------------
    NXS <- ncol( XS)
    if(is.null(colnames(XS)))
        colnames(XS) <- rep("XS", NXS)
    NXO <- ncol( XO)
    if(is.null(colnames(XO)))
        colnames(XO) <- rep("XO", NXO)
    Nparam <- NXS + NXO + 2
                                        # Total # of parameters
    nObs <- length( YS)
    NO <- length( YS[YS > 0])
   nParam <- NXS + NXO + 2
   ## parameter indices
   iBetaS <- 1:NXS
   iBetaO <- seq(tail(iBetaS, 1)+1, length=NXO)
   iSigma <- tail(iBetaO, 1) + 1
   iRho <- tail(iSigma, 1) + 1
   ## output, if asked for it
   if( print.level > 0) {
      cat( "Yo observed:", NO, "times; not observed:", nObs - NO, "times\n")
      cat( "Initial values:\n")
      cat("Selection\n")
      print(start[iBetaS])
      cat("Outcome\n")
      print(start[iBetaO])
      cat("sigma1\n")
      print(start[iSigma])
      cat("rho1\n")
      print(start[iRho])
   }
   ## pre-calculate a few values:
   XS0 <- XS[YS==0,,drop=FALSE]
   XS1 <- XS[YS==1,,drop=FALSE]
   YO1 <- YO[YS==1]
   XO1 <- XO[YS==1,,drop=FALSE]
   N0 <- sum(YS==0)
   N1 <- sum(YS==1)
   ## estimate
   result <- maxLik(loglik, grad=gradlik, hess=hesslik, start=start,
                     print.level=print.level,
                    method=maxMethod, ...)
   result$tobitType <- 2
   result$method <- "ml"
   class( result ) <- c( "selection", class( result ) )
   return( result )
}
