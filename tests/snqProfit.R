library( micEcon )
data( germanFarms )
options( digits = 6 )

germanFarms$qOutput   <- germanFarms$vOutput / germanFarms$pOutput
germanFarms$qVarInput <- -germanFarms$vVarInput / germanFarms$pVarInput
germanFarms$qLabor    <- -germanFarms$qLabor
germanFarms$time      <- c( 0:19 )

pNamesT <- c( "pOutput", "pVarInput", "pLabor" )
qNamesT <- c( "qOutput", "qVarInput", "qLabor" )
fNamesT <- c( "land", "time" )

estResult <- snqProfitEst( pNamesT, qNamesT, "land", data = germanFarms )
print( estResult )
class( estResult ) <- NULL
estResult$est <- summary( estResult$est )
print( estResult )

################ without fix inputs ##############################
estResult <- snqProfitEst( pNamesT, qNamesT, NULL, data = germanFarms )
print( estResult )
class( estResult ) <- NULL
estResult$est <- summary( estResult$est )
print( estResult )

estResultCalc <- snqProfitCalc( pNamesT, NULL, estResult$data,
   estResult$weights, estResult$scalingFactors, estResult$coef )
print( estResultCalc )
if( max( abs( estResultCalc - estResult$fitted ) ) > 1e-5 ) {
   stop( "values from snqProfitCalc are not equal to fitted values" )
}

estResultEla <- snqProfitEla( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$data[ 20, qNamesT ],
   estResult$weights, estResult$scalingFactors )
print( estResultEla )

estResultHessianDeriv <- snqProfitHessianDeriv( estResult$pMean,
   estResult$weights, nFix = 2 )
print( estResultHessianDeriv )

estResultHessian <- snqProfitHessian( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$weights,
   estResult$scalingFactors )
print( estResultHessian )

########### with fix inputs, form = 0 ########################
estResult <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms )
print( estResult )
class( estResult ) <- NULL
estResult$est <- summary( estResult$est )
print( estResult )

estResultCalc <- snqProfitCalc( pNamesT, fNamesT, estResult$data,
   estResult$weights, estResult$scalingFactors, estResult$coef )
print( estResultCalc )
if( max( abs( estResultCalc - estResult$fitted ) ) > 1e-5 ) {
   stop( "values from snqProfitCalc are not equal to fitted values" )
}

estResultEla <- snqProfitEla( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$data[ 20, qNamesT ],
   estResult$weights, estResult$scalingFactors )
print( estResultEla )

estResultHessianDeriv <- snqProfitHessianDeriv( estResult$pMean,
   estResult$weights, nFix = 2 )
print( estResultHessianDeriv )

estResultHessian <- snqProfitHessian( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$weights,
   estResult$scalingFactors )
print( estResultHessian )

estResultShadowprices <- snqProfitShadowPrices( pNamesT, fNamesT, estResult )
print( estResultShadowprices )

####################################################
estResult <- snqProfitEst( pNamesT, qNamesT, fNamesT, data=germanFarms, form = 1 )
print( estResult )
class( estResult ) <- NULL
estResult$est <- summary( estResult$est )
print( estResult )

estResultCalc <- snqProfitCalc( pNamesT, fNamesT, estResult$data,
   estResult$weights, estResult$scalingFactors, estResult$coef, form = 1 )
print( estResultCalc )
if( max( abs( estResultCalc - estResult$fitted ) ) > 1e-5 ) {
   stop( "values from snqProfitCalc are not equal to fitted values" )
}

estResultEla <- snqProfitEla( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$data[ 20, qNamesT ],
   estResult$weights, estResult$scalingFactors )
print( estResultEla )

estResultHessianDeriv <- snqProfitHessianDeriv( estResult$pMean,
   estResult$weights, nFix = 2, form = 1 )
print( estResultHessianDeriv )

estResultHessian <- snqProfitHessian( estResult$coef$beta,
   estResult$data[ 20, pNamesT ], estResult$weights,
   estResult$scalingFactors )
print( estResultHessian )

estResultShadowprices <- snqProfitShadowPrices( pNamesT, fNamesT, estResult )
print( estResultShadowprices )
