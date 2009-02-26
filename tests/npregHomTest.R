# load micEcon package
library( micEcon )

# load data
data( germanFarms )
# output quantity
germanFarms$qOutput <- germanFarms$vOutput / germanFarms$pOutput
# quantity of variable inputs
germanFarms$qVarInput <- germanFarms$vVarInput / germanFarms$pVarInput
# a time trend to account for technical progress:
germanFarms$time <- c(1:20)

# weights to impose normalize prices
weights <- c(
   pOutput = mean( germanFarms$qOutput ),
   pVarInput = mean( germanFarms$qVarInput ),
   pLabor = mean( germanFarms$qLabor ) )
weights <- weights / sum( weights )

# estimation (restricted gradients)
npseed( 123 )
estResult <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights )
print( estResult )
all.equal( estResult$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResult$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResult$grad[ , "pLabor" ] * germanFarms$pLabor )
estEla <- elas( estResult )
print( estEla )
all.equal( estEla[[ "pOutput" ]] + estEla[[ "pVarInput" ]],
   - estEla[[ "pLabor" ]] )
estElaObs <- elas( estResult, yObs = TRUE )
print( estElaObs )
all.equal( estElaObs[[ "pOutput" ]] + estElaObs[[ "pVarInput" ]],
   - estElaObs[[ "pLabor" ]] )
# different normalized variable omitted
npseed( 123 )
estResult2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ] )
print( estResult2 )
all.equal( estResult2$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResult2$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResult2$grad[ , "pLabor" ] * germanFarms$pLabor )
estEla2 <- elas( estResult2 )
print( estEla2 )
all.equal( estEla2[[ "pOutput" ]] + estEla2[[ "pVarInput" ]],
   - estEla2[[ "pLabor" ]] )
estEla2Obs <- elas( estResult2, yObs = TRUE )
print( estEla2Obs )
all.equal( estEla2Obs[[ "pOutput" ]] + estEla2Obs[[ "pVarInput" ]],
   - estEla2Obs[[ "pLabor" ]] )

# estimation (gradients not restricted)
npseed( 123 )
estResultAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE )
print( estResultAll )
all.equal( estResultAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultAll$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaAll <- elas( estResultAll )
print( estElaAll )
all.equal( estElaAll[[ "pOutput" ]] + estElaAll[[ "pVarInput" ]],
   - estElaAll[[ "pLabor" ]] )
estElaAllObs <- elas( estResultAll, yObs = TRUE )
print( estElaAllObs )
all.equal( estElaAllObs[[ "pOutput" ]] + estElaAllObs[[ "pVarInput" ]],
   - estElaAllObs[[ "pLabor" ]] )
# different order of weights
npseed( 123 )
estResultAll2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ],
   restrictGrad = FALSE )
all.equal( estResultAll$grad, estResultAll2$grad, tolerance = 1e-6 )
all.equal( elas( estResultAll ), elas( estResultAll2 ), tolerance = 1e-6 )
all.equal( elas( estResultAll, yObs = TRUE ),
   elas( estResultAll2, yObs = TRUE ), tolerance = 1e-6 )


# estimation with Epanechnikov kernel (restricted gradients)
npseed( 123 )
estResultEpa <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, ckertype="epanechnikov" )
print( estResultEpa )
all.equal( estResultEpa$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultEpa$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultEpa$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaEpa <- elas( estResultEpa )
print( estElaEpa )
all.equal( estElaEpa[[ "pOutput" ]] + estElaEpa[[ "pVarInput" ]],
   - estElaEpa[[ "pLabor" ]] )
estElaEpaObs <- elas( estResultEpa, yObs = TRUE )
print( estElaEpaObs )
all.equal( estElaEpaObs[[ "pOutput" ]] + estElaEpaObs[[ "pVarInput" ]],
   - estElaEpaObs[[ "pLabor" ]] )

# estimation with Epanechnikov kernel (gradients not restricted)
npseed( 123 )
estResultEpaAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE,
   ckertype="epanechnikov" )
print( estResultEpaAll )
all.equal( estResultEpaAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultEpaAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultEpaAll$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaEpaAll <- elas( estResultEpaAll )
print( estElaEpaAll )
all.equal( estElaEpaAll[[ "pOutput" ]] + estElaEpaAll[[ "pVarInput" ]],
   - estElaEpaAll[[ "pLabor" ]] )
estElaEpaAllObs <- elas( estResultEpaAll, yObs = TRUE )
print( estElaEpaAllObs )
all.equal( estElaEpaAllObs[[ "pOutput" ]] + estElaEpaAllObs[[ "pVarInput" ]],
   - estElaEpaAllObs[[ "pLabor" ]] )
# different order of weights
npseed( 123 )
estResultEpaAll2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ],
   restrictGrad = FALSE, ckertype="epanechnikov" )
all.equal( estResultEpaAll$grad, estResultEpaAll2$grad, tolerance = 1e-6 )
all.equal( elas( estResultEpaAll ), elas( estResultEpaAll2 ), tolerance = 1e-6 )
all.equal( elas( estResultEpaAll, yObs = TRUE ),
   elas( estResultEpaAll2, yObs = TRUE ), tolerance = 1e-6 )


# estimation with manual bandwidth selection (restricted gradients)
estResultMan <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, bws = rep( 1, 3 ),
   bwscaling = TRUE )
print( estResultMan )
all.equal( estResultMan$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultMan$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultMan$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaMan <- elas( estResultMan )
print( estElaMan )
all.equal( estElaMan[[ "pOutput" ]] + estElaMan[[ "pVarInput" ]],
   - estElaMan[[ "pLabor" ]] )
estElaManObs <- elas( estResultMan, yObs = TRUE )
print( estElaManObs )
all.equal( estElaManObs[[ "pOutput" ]] + estElaManObs[[ "pVarInput" ]],
   - estElaManObs[[ "pLabor" ]] )
# different normalized variable omitted
estResultMan2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ],
   bws = rep( 1, 3 ), bwscaling = TRUE )
print( estResultMan2 )
all.equal( estResultMan2$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultMan2$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultMan2$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaMan2 <- elas( estResultMan2 )
print( estElaMan2 )
all.equal( estElaMan2[[ "pOutput" ]] + estElaMan2[[ "pVarInput" ]],
   - estElaMan2[[ "pLabor" ]] )
estElaMan2Obs <- elas( estResultMan2, yObs = TRUE )
print( estElaMan2Obs )
all.equal( estElaMan2Obs[[ "pOutput" ]] + estElaMan2Obs[[ "pVarInput" ]],
   - estElaMan2Obs[[ "pLabor" ]] )

# estimation with manual bandwidth selection (gradients not restricted)
estResultManAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE,
   bws = rep( 1, 4 ), bwscaling = TRUE )
print( estResultManAll )
all.equal( estResultManAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultManAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultManAll$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaManAll <- elas( estResultManAll )
print( estElaManAll )
all.equal( estElaManAll[[ "pOutput" ]] + estElaManAll[[ "pVarInput" ]],
   - estElaManAll[[ "pLabor" ]] )
estElaManAllObs <- elas( estResultManAll, yObs = TRUE )
print( estElaManAllObs )
all.equal( estElaManAllObs[[ "pOutput" ]] + estElaManAllObs[[ "pVarInput" ]],
   - estElaManAllObs[[ "pLabor" ]] )
# different order of weights
estResultManAll2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ],
   restrictGrad = FALSE, bws = rep( 1, 4 ), bwscaling = TRUE )
all.equal( estResultManAll$grad, estResultManAll2$grad )
all.equal( elas( estResultManAll ), elas( estResultManAll2 ) )
all.equal( elas( estResultManAll, yObs = TRUE ),
   elas( estResultManAll2, yObs = TRUE ) )


# local-linear estimation (restricted gradients)
npseed( 123 )
estResultLl <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, regtype = "ll" )
print( estResultLl )
all.equal( estResultLl$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultLl$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultLl$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaLl <- elas( estResultLl )
print( estElaLl )
all.equal( estElaLl[[ "pOutput" ]] + estElaLl[[ "pVarInput" ]],
   - estElaLl[[ "pLabor" ]] )
estElaLlObs <- elas( estResultLl, yObs = TRUE )
print( estElaLlObs )
all.equal( estElaLlObs[[ "pOutput" ]] + estElaLlObs[[ "pVarInput" ]],
   - estElaLlObs[[ "pLabor" ]] )

# local-linear estimation (gradients not restricted)
npseed( 123 )
estResultLlAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE,
   regtype = "ll" )
print( estResultLlAll )
all.equal( estResultLlAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultLlAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultLlAll$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaLlAll <- elas( estResultLlAll )
print( estElaLlAll )
all.equal( estElaLlAll[[ "pOutput" ]] + estElaLlAll[[ "pVarInput" ]],
   - estElaLlAll[[ "pLabor" ]] )
estElaLlAllObs <- elas( estResultLlAll, yObs = TRUE )
print( estElaLlAllObs )
all.equal( estElaLlAllObs[[ "pOutput" ]] + estElaLlAllObs[[ "pVarInput" ]],
   - estElaLlAllObs[[ "pLabor" ]] )


# local-linear estimation with Epanechnikov kernel (restricted gradients)
npseed( 123 )
estResultLlEpa <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, ckertype="epanechnikov",
   regtype = "ll" )
print( estResultLlEpa )
all.equal( estResultLlEpa$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultLlEpa$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultLlEpa$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaLlEpa <- elas( estResultLlEpa )
print( estElaLlEpa )
all.equal( estElaLlEpa[[ "pOutput" ]] + estElaLlEpa[[ "pVarInput" ]],
   - estElaLlEpa[[ "pLabor" ]] )
estElaLlEpaObs <- elas( estResultLlEpa, yObs = TRUE )
print( estElaLlEpaObs )
all.equal( estElaLlEpaObs[[ "pOutput" ]] + estElaLlEpaObs[[ "pVarInput" ]],
   - estElaLlEpaObs[[ "pLabor" ]] )

# local linear estimation with Epanechnikov kernel (gradients not restricted)
npseed( 123 )
estResultLlEpaAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE,
   ckertype = "epanechnikov", regtype = "ll" )
print( estResultLlEpaAll )
all.equal( estResultLlEpaAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultLlEpaAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultLlEpaAll$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaLlEpaAll <- elas( estResultLlEpaAll )
print( estElaLlEpaAll )
all.equal( estElaLlEpaAll[[ "pOutput" ]] + estElaLlEpaAll[[ "pVarInput" ]],
   - estElaLlEpaAll[[ "pLabor" ]] )
estElaLlEpaAllObs <- elas( estResultLlEpaAll, yObs = TRUE )
print( estElaLlEpaAllObs )
all.equal( estElaLlEpaAllObs[[ "pOutput" ]] + estElaLlEpaAllObs[[ "pVarInput" ]],
   - estElaLlEpaAllObs[[ "pLabor" ]] )


# local-linear estimation with manual bandwidth selection (restricted gradients)
estResultLlMan <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, bws = rep( 1, 3 ),
   bwscaling = TRUE, regtype = "ll" )
print( estResultLlMan )
all.equal( estResultLlMan$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultLlMan$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultLlMan$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaLlMan <- elas( estResultLlMan )
print( estElaLlMan )
all.equal( estElaLlMan[[ "pOutput" ]] + estElaLlMan[[ "pVarInput" ]],
   - estElaLlMan[[ "pLabor" ]] )
estElaLlManObs <- elas( estResultLlMan, yObs = TRUE )
print( estElaLlManObs )
all.equal( estElaLlManObs[[ "pOutput" ]] + estElaLlManObs[[ "pVarInput" ]],
   - estElaLlManObs[[ "pLabor" ]] )

# local-linear estimation with manual bandwidth selection (gradients not restricted)
estResultLlManAll <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights, restrictGrad = FALSE,
   bws = rep( 1, 4 ), bwscaling = TRUE, regtype = "ll" )
print( estResultLlManAll )
all.equal( estResultLlManAll$grad[ , "pOutput" ] * germanFarms$pOutput +
   estResultLlManAll$grad[ , "pVarInput" ] * germanFarms$pVarInput,
   - estResultLlManAll$grad[ , "pLabor" ] * germanFarms$pLabor )
estElaLlManAll <- elas( estResultLlManAll )
print( estElaLlManAll )
all.equal( estElaLlManAll[[ "pOutput" ]] + estElaLlManAll[[ "pVarInput" ]],
   - estElaLlManAll[[ "pLabor" ]] )
estElaLlManAllObs <- elas( estResultLlManAll, yObs = TRUE )
print( estElaLlManAllObs )
all.equal( estElaLlManAllObs[[ "pOutput" ]] + estElaLlManAllObs[[ "pVarInput" ]],
   - estElaLlManAllObs[[ "pLabor" ]] )
# different order of weights
estResultLlManAll2 <- npregHom( "qVarInput",
   xNames = c( "pOutput", "pVarInput", "pLabor", "land" ),
   data = germanFarms, homWeights = weights[ c( 3, 2, 1 ) ], restrictGrad = FALSE,
   bws = rep( 1, 4 ), bwscaling = TRUE, regtype = "ll" )
all.equal( estResultLlManAll$grad, estResultLlManAll2$grad )
all.equal( estResultLlManAll$est, estResultLlManAll2$est, check.attributes = FALSE )
