snqProfitImposeConvexity <- function( estResult, rankReduction = 0,
   start = 10, ... ) {

   if( class( estResult ) != "snqProfitEst" ) {
      stop( "argument 'estResult' must be of class 'snqProfitEst'." )
   }
   pNames  <- names( estResult$pMeans )
   qNames  <- names( estResult$qMeans )
   fNames  <- names( estResult$fMeans )
   nNetput <- length( pNames )
   nFix    <- length( fNames )
   result  <- list()

   ## preparations for minimum distance estimation
   uVecliHessian <- vecli( estResult$hessian[ 1:( nNetput - 1 ),
      1:( nNetput - 1 ) ] )
      # vector of linear indep. values of unconstr. Hessian
   hessianDeriv <- snqProfitHessianDeriv( estResult$pMeans, estResult$weights,
      nFix = nFix, form = estResult$form )
      # derivatives of the Hessian with respect to the coefficients

   ## distance function between the unconstrained and constrained Hessian
   hessianDistance <- function( cholVec ) {
         # cholVec = vector of values of triangular Cholesky Matrix
      cholMat <- triang( cholVec, nNetput - 1 )     # triangular Cholesky Matrix
      cVecliHessian <- vecli( t(cholMat) %*% cholMat )
         # vector of linear independent values of constrained Hessian
      result <- t( uVecliHessian - cVecliHessian ) %*% solve( hessianDeriv %*%
         estResult$coef$liCoefCov %*% t( hessianDeriv ) ) %*%
         ( uVecliHessian - cVecliHessian )
      return( result )
   }

   ## non-linear minimization of the distance function
   nCholValues <- nNetput * ( nNetput - 1 ) / 2 -
      rankReduction * ( rankReduction + 1 ) / 2
      # number of non-zero values in the Cholesky matrix
   cholVec <- array( start, c( nCholValues ) ) # starting values
   mindist <- optim( cholVec, hessianDistance, ... )
   if( mindist$convergence != 0 ) {
      if( mindist$convergence == 1 ) {
         stop( "non-linear minimization with optim(): iteration limit exceeded." )
      } else {
         stop( "non-linear minimization: 'optim' did not converge." )
      }
   }

   ## asymptotic least squares (ALS)
   ## (all coefficients may adjust to best fit of the model)
   cholMat <- triang( mindist$par, nNetput - 1 )
      # triangular Cholesky Matrix
   cVecliHessian <- vecli( t(cholMat) %*% cholMat )
      # vector of linear independent values of constrained Hessian
   coef <- estResult$coef$liCoef + estResult$coef$liCoefCov %*%
      t( hessianDeriv ) %*% solve( hessianDeriv %*% estResult$coef$liCoefCov %*%
      t( hessianDeriv ) ) %*% ( cVecliHessian - uVecliHessian )
      # vector of li indep. constrained coefficients

   ## results of constrained model
   result$pMeans <- estResult$pMeans
   result$qMeans <- estResult$qMeans
   result$fMeans <- estResult$fMeans
   result$mindist <- mindist
   result$coef <- snqProfitCoef( coef, nNetput, nFix, form = estResult$form,
      qNames = names( estResult$qMeans ), pNames = names( estResult$pMeans ),
      fNames = names( estResult$fMeans ) )
      # constrained coefficients
   result$fitted <- snqProfitCalc( pNames, fNames, data = estResult$estData,
      weights = estResult$weights, coef = result$coef, form = estResult$form )
   result$residuals <- data.frame( nr = c( 1:nrow( estResult$estData ) ) )
   for( i in 1:nNetput ) {
      result$residuals[[ qNames[ i ] ]] <- estResult$estData[[ qNames[ i ] ]] -
         result$fitted[ , i ]
   }
   if( !( "nr" %in% qNames ) ) {
      result$residuals[[ "nr" ]] <- NULL
   }
   result$r2 <- array( NA, c( nNetput ) )
   for( i in 1:nNetput ) {
      result$r2[ i ] <- rSquared( estResult$estData[[ qNames[ i ] ]],
         result$residuals[[ qNames[ i ] ]] )
   }
   names( result$r2 ) <- names( estResult$qMeans )

   result$hessian <- snqProfitHessian( result$coef$beta, estResult$pMean,
      estResult$weights ) # constrained Hessian matrix
   result$ela <- snqProfitEla( result$coef$beta, estResult$pMean, estResult$qMean,
      estResult$weights ) # elasticities of constrained model
   result$estData <- estResult$estData
   result$weights <- estResult$weights
   result$normPrice <- estResult$normPrice
   result$convexity <- TRUE
   result$form      <- estResult$form

   class( result ) <- "snqProfitImposeConvexity"
   return( result )
}