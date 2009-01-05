library( micEcon )
library( plm )

## preparing data
data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

## estimate a quadratic production function
estResult <- quadFuncEst( "qOutput",
   c( "qLabor", "land", "qVarInput", "time" ), germanFarms )
coef( estResult )
print( estResult )

## compute fitted values
fitted <- quadFuncCalc( c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ) )
all.equal( fitted, estResult$fitted )

## compute the marginal products of the inputs
margProducts <- quadFuncDeriv(
   c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, coef( estResult ), vcov( estResult ) )
print( margProducts )

## estimate a quadratic production function with a shifter
estResultShifter <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = "time", data = germanFarms )
coef( estResultShifter )
print( estResultShifter )

## estimate a quadratic production function with 2 shifters
germanFarms$timeSq <- germanFarms$time^2
estResultShifter2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput" ),
   shifterNames = c( "time", "timeSq" ), data = germanFarms )
coef( estResultShifter2 )
print( estResultShifter2 )

## estimate a linear functions with quadFuncEst
estResultLinear <- quadFuncEst( yName = "qOutput", xNames = NULL,
   shifterNames = c( "time", "qLabor", "land", "qVarInput" ),
   data = germanFarms )
coef( estResultLinear )
print( estResultLinear )

## estimate with further argument passed to lm()
estResult2 <- quadFuncEst( yName = "qOutput",
   xNames = c( "qLabor", "land", "qVarInput", "time" ),
   germanFarms, x = TRUE, y = TRUE )
coef( estResult2 )
print( estResult2 )

## panel data
data( "GrunfeldGreene", package = "systemfit" )
ggData <- plm.data( GrunfeldGreene, c( "firm", "year" ) )
# fixed effects
ggResult <- quadFuncEst( "invest", c( "value", "capital" ), ggData )
coef( ggResult )
print( ggResult )
# random effects
ggResultRan <- quadFuncEst( "invest", c( "value", "capital" ), ggData,
   model = "random", random.method = "amemiya" )
coef( ggResultRan )
print( ggResultRan )
