library( micEcon )

data( germanFarms )
# output quantity:
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

# estimate a Cobb-Douglas production function
estResult <- translogEst( "qOutput", c( "qLabor", "qVarInput", "land", "time" ),
   germanFarms, linear = TRUE )

# calculate fitted values
fitted <- cobbDouglasCalc( c( "qLabor", "qVarInput", "land", "time" ),
   data = germanFarms, coef = coef( estResult )[ 1:5 ] )
print( fitted )
all.equal( fitted, estResult$fitted )

# calculate fitted values using logged independent variables
germanFarms$lQLabor    <- log( germanFarms$qLabor )
germanFarms$lLand      <- log( germanFarms$land )
germanFarms$lQVarInput <- log( germanFarms$qVarInput )
germanFarms$lTime      <- log( germanFarms$time )
fittedLogged <- cobbDouglasCalc( c( "lQLabor", "lQVarInput", "lLand", "lTime" ),
   data = germanFarms, coef = coef( estResult )[ 1:5 ], dataLogged = TRUE )
all.equal( fitted, exp( fittedLogged ) )

# coefficients not named
coefNoNames <- coef( estResult )[ 1:5 ]
names( coefNoNames ) <- NULL
fittedNoNames <- cobbDouglasCalc( c( "qLabor", "qVarInput", "land", "time" ),
   data = germanFarms, coef = coefNoNames )
all.equal( fitted, fittedNoNames )

# coefficients in a different order
coefDiffOrder <- coef( estResult )[ c( 3, 5, 1, 2, 4 ) ]
fittedDiffOrder <- cobbDouglasCalc( c( "qLabor", "qVarInput", "land", "time" ),
   data = germanFarms, coef = coefDiffOrder )
all.equal( fitted, fittedDiffOrder )

# calculate optimal quantities of variable inputs
xCoef <- coef( estResult )[ 1:3 ]
zCoef <- coef( estResult )[ 4:5 ]
names( zCoef ) <- c( "d_1", "d_2" )
optInput <- cobbDouglasOpt( pyName = "pOutput",
   pxNames = c( "pLabor", "pVarInput" ), coef = xCoef,
   data = germanFarms, xNames = c( "qLabor", "qVarInput" ),
   zNames = c( "land", "time" ), zCoef = zCoef )
print( optInput )

# determine optimal quantities of variable inputs using optim()
objFun <- function( xVal, obs = 1 ) {
   tmpData <- germanFarms
   tmpData$qLabor[ obs ] <- xVal[ 1 ]
   tmpData$qVarInput[ obs ] <- xVal[ 2 ]
   outp <- translogCalc( c( "qLabor", "qVarInput", "land", "time" ),
      data = tmpData, coef = coef( estResult ) )
   profit <- germanFarms$pOutput[ obs ] * outp[ obs ] -
      germanFarms$pLabor[ obs ] * xVal[ 1 ] -
      germanFarms$pVarInput[ obs ] * xVal[ 2 ]
   return( profit )
}
optInputNum <- data.frame( qLabor = rep( NA, nrow( germanFarms ) ),
   qVarInput = rep( NA, nrow( germanFarms ) ) )
for( obs in 1:nrow( germanFarms ) ) {
   optResult <- optim(
      c( germanFarms$qLabor[ obs ], germanFarms$qVarInput[ obs ] ),
      objFun, method = "L-BFGS-B", lower = 1e-10,
      control = list( fnscale = -1 ), obs = obs )
   optInputNum[ obs, ] <- optResult$par
}
all.equal( optInput, optInputNum, check.attributes = FALSE, tolerance = 1e-5 )
