library( micEcon )

## preparing data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

## testing translogEst
# estimate a translog production function
estResult <- translogEst( "qOutput", c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms )

print( estResult )

## testing translogCalc
fitted <- translogCalc( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef )

all.equal( fitted, estResult$fitted )

## testing translogDeriv
margProducts <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, estResult$coefCov )
print( margProducts )

margProductsObs <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, estResult$coefCov, yName = "qOutput" )
print( margProductsObs )

germanFarms$fitted <- fitted
margProductsFitted <- translogDeriv( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, estResult$coefCov, yName = "fitted" )
all.equal( margProducts, margProductsFitted )


## testing translogHessian
# compute the Hessian matrices
hessians <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef )
print( hessians )

hessiansObs <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, yName = "qOutput" )
print( hessiansObs )

germanFarms$fitted <- fitted
hessiansFitted <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, yName = "fitted" )
all.equal( hessians, hessiansFitted )

# compute the bordered Hessian matrices
borderedHessians <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, bordered = TRUE )
print( borderedHessians )

borderedHessiansObs <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, yName = "qOutput", bordered = TRUE )
print( borderedHessiansObs )

germanFarms$fitted <- fitted
borderedHessiansFitted <- translogHessian( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, estResult$coef, yName = "fitted", bordered = TRUE )
all.equal( borderedHessians, borderedHessiansFitted )

## testing translogCheckMono
test <- translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ) )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), strict = TRUE )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), increasing = FALSE )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckMono( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), increasing = FALSE, strict = TRUE )
summary( test )
class( test ) <- NULL
print( test )


## testing translogCheckCurvature
test <- translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ) )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), quasi = TRUE )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), convexity = FALSE )
summary( test )
class( test ) <- NULL
print( test )

test <- translogCheckCurvature( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), convexity = FALSE, quasi = TRUE )
summary( test )
class( test ) <- NULL
print( test )
