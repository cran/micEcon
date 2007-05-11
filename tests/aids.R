library( micEcon )
data( Blanciforti86 )
options( digits = 3 )

set <- !is.na( Blanciforti86$pFood1 )
setWo1 <- set & rownames( Blanciforti86 ) != 1947
pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
pMeans <- colMeans( Blanciforti86[ set, pNames ] )

cat( paste( "\nRepeating the demand analysis of Blanciforti, Green",
   "& King (1986)\n" ) )
estResultLA <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "LA:SL" )
print( estResultLA )
print( summary( estResultLA ) )
print( elas( estResultLA, method = "Ch", quantNames = wNames ) )
# imposing restrictions via TX
estResultLATX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "LA:SL",
   TX = TRUE )
print( estResultLATX )
print( summary( estResultLATX ) )
print( elas( estResultLATX, method = "Ch", quantNames = wNames ) )
estResultLATX$call <- NULL
estResultLATX$est$bt <- NULL
estResultLATX$est$btcov <- NULL
estResultLATX$est$x <- NULL
estResultLATX$est$TX <- NULL
estResultLATX$est$q.restr <- NULL
estResultLA$call <- NULL
estResultLA$est$x <- NULL
estResultLA$est$R.restr <- NULL
estResultLA$est$q.restr <- NULL
print( all.equal( estResultLA, estResultLATX ) )

## only homogeneity (no symmetry imposed)
estResultLAhom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL" )
print( estResultLAhom )
print( summary( estResultLAhom ) )
print( elas( estResultLAhom, method = "Ch", quantNames = wNames ) )
# imposing restrictions via TX
estResultLAhomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL",
   TX = TRUE )
print( estResultLAhomTX )
print( summary( estResultLAhomTX ) )
print( elas( estResultLAhomTX, method = "Ch", quantNames = wNames ) )
estResultLAhomTX$call <- NULL
estResultLAhomTX$est$bt <- NULL
estResultLAhomTX$est$btcov <- NULL
estResultLAhomTX$est$x <- NULL
estResultLAhomTX$est$TX <- NULL
estResultLAhomTX$est$q.restr <- NULL
estResultLAhom$call <- NULL
estResultLAhom$est$x <- NULL
estResultLAhom$est$R.restr <- NULL
estResultLAhom$est$q.restr <- NULL
print( all.equal( estResultLAhom, estResultLAhomTX ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultLAunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL" )
print( estResultLAunr )
print( summary( estResultLAunr ) )
print( elas( estResultLAunr, method = "Ch", quantNames = wNames ) )
# imposing restrictions via TX
estResultLAunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL",
   TX = TRUE )
print( estResultLAunrTX )
print( summary( estResultLAunrTX ) )
print( elas( estResultLAunrTX, method = "Ch", quantNames = wNames ) )
estResultLAunrTX$call <- NULL
estResultLAunrTX$est$bt <- NULL
estResultLAunrTX$est$btcov <- NULL
estResultLAunrTX$est$x <- NULL
estResultLAunrTX$est$TX <- NULL
estResultLAunrTX$est$q.restr <- NULL
estResultLAunr$call <- NULL
estResultLAunr$est$x <- NULL
estResultLAunr$est$R.restr <- NULL
estResultLAunr$est$q.restr <- NULL
print( all.equal( estResultLAunr, estResultLAunrTX ) )

#####################################################
cat( paste( "\nRepeating the evaluation of different elasticity formulas",
   "of Green & Alston (1990): iterated AIDS\n" ) )
estResultAIDS <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], method = "IL:L" )
print( estResultAIDS )
print( summary( estResultAIDS ) )
print( elas( estResultAIDS, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via TX
estResultAIDSTX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], method = "IL:L", TX = TRUE )
print( estResultAIDSTX )
print( summary( estResultAIDSTX ) )
print( elas( estResultAIDSTX, method = "AIDS", quantNames = wNames ) )
estResultAIDSTX$call <- NULL
estResultAIDSTX$est$bt <- NULL
estResultAIDSTX$est$btcov <- NULL
estResultAIDSTX$est$x <- NULL
estResultAIDSTX$est$TX <- NULL
estResultAIDSTX$est$q.restr <- NULL
estResultAIDS$call <- NULL
estResultAIDS$est$x <- NULL
estResultAIDS$est$R.restr <- NULL
estResultAIDS$est$q.restr <- NULL
print( all.equal( estResultAIDS, estResultAIDSTX ) )

## only homogeneity (no symmetry imposed)
estResultAIDShom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL:L" )
print( estResultAIDShom )
print( summary( estResultAIDShom ) )
print( elas( estResultAIDShom, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via TX
estResultAIDShomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL:L", TX = TRUE )
print( estResultAIDShomTX )
print( summary( estResultAIDShomTX ) )
print( elas( estResultAIDShomTX, method = "AIDS", quantNames = wNames ) )
estResultAIDShomTX$call <- NULL
estResultAIDShomTX$est$bt <- NULL
estResultAIDShomTX$est$btcov <- NULL
estResultAIDShomTX$est$x <- NULL
estResultAIDShomTX$est$TX <- NULL
estResultAIDShomTX$est$q.restr <- NULL
estResultAIDShom$call <- NULL
estResultAIDShom$est$x <- NULL
estResultAIDShom$est$R.restr <- NULL
estResultAIDShom$est$q.restr <- NULL
print( all.equal( estResultAIDShom, estResultAIDShomTX ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultAIDSunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL:L" )
print( estResultAIDSunr )
print( summary( estResultAIDSunr ) )
print( elas( estResultAIDSunr, method = "AIDS", quantNames = wNames ) )
# imposing restrictions via TX
estResultAIDSunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], method = "IL:L", TX = TRUE )
print( estResultAIDSunrTX )
print( summary( estResultAIDSunrTX ) )
print( elas( estResultAIDSunrTX, method = "AIDS", quantNames = wNames ) )
estResultAIDSunrTX$call <- NULL
estResultAIDSunrTX$est$bt <- NULL
estResultAIDSunrTX$est$btcov <- NULL
estResultAIDSunrTX$est$x <- NULL
estResultAIDSunrTX$est$TX <- NULL
estResultAIDSunrTX$est$q.restr <- NULL
estResultAIDSunr$call <- NULL
estResultAIDSunr$est$x <- NULL
estResultAIDSunr$est$R.restr <- NULL
estResultAIDSunr$est$q.restr <- NULL
print( all.equal( estResultAIDSunr, estResultAIDSunrTX ) )

## with NAs
estResultLaSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   method = "LA:S" )
print( estResultLaSNa )
print( summary( estResultLaSNa ) )
print( elas( estResultLaSNa, method = "AIDS", quantNames = wNames ) )

estResultLaSlNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   method = "LA:SL" )
print( estResultLaSlNa )
print( summary( estResultLaSlNa ) )
print( elas( estResultLaSlNa, method = "AIDS", quantNames = wNames ) )

estResultLaLNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86,
   method = "LA:L" )
print( estResultLaLNa )
print( summary( estResultLaLNa ) )
print( elas( estResultLaLNa, method = "AIDS", quantNames = wNames ) )

estResultAIDSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, method = "IL:L" )
print( estResultAIDSNa )
print( summary( estResultAIDSNa ) )
print( elas( estResultAIDSNa, method = "AIDS", quantNames = wNames ) )


########## Elasticities ###############
cat( "\nAIDS: Elasticities\n" )
ela <- aidsElas( estResultAIDS$coef, wMeans, pMeans, method = "AIDS",
   coefVcov = estResultAIDS$coef$allcov, df = estResultAIDS$est$df )
print( ela )
print( summary( ela ) )
elaTX <- aidsElas( estResultAIDSTX$coef, wMeans, pMeans, method = "AIDS",
   coefVcov = estResultAIDSTX$coef$allcov, df = estResultAIDSTX$est$df )
print( elaTX )
print( summary( elaTX ) )
print( all.equal( ela, elaTX ) )

print( elas( estResultAIDS ) )
print( summary( elas( estResultAIDS ) ) )

print( elas( estResultAIDSTX ) )
print( summary( elas( estResultAIDSTX ) ) )


cat( "\nLA: Elasticity formula of non-linear AIDS\n" )
ela <- aidsElas( estResultLA$coef, wMeans, pMeans, method = "AIDS",
   coefVcov = estResultLA$coef$allcov, df = estResultLA$est$df )
print( ela )
print( summary( ela ) )
elaTX <- aidsElas( estResultLATX$coef, wMeans, pMeans, method = "AIDS",
   coefVcov = estResultLATX$coef$allcov, df = estResultLATX$est$df )
print( elaTX )
print( summary( elaTX ) )
print( all.equal( ela, elaTX ) )

print( elas( estResultLA, method = "AIDS" ) )
print( summary( elas( estResultLA, method = "AIDS" ) ) )

print( elas( estResultLATX, method = "AIDS" ) )
print( summary( elas( estResultLATX, method = "AIDS" ) ) )


cat( "\n********** Elasticities ***************" )
cat( "\nLA: Elasticity formula of Goddard or Chalfant\n" )
ela <- aidsElas( estResultLA$coef, wMeans, method = "Go",
   coefVcov = estResultLA$coef$allcov, df = estResultLA$est$df )
print( ela )
print( summary( ela ) )
ela <- aidsElas( estResultLA$coef, wMeans, method = "Ch",
   coefVcov = estResultLA$coef$allcov, df = estResultLA$est$df )
print( ela )
print( summary( ela ) )

print( elas( estResultLA, method = "Go" ) )
print( summary( elas( estResultLA ) ) )

print( elas( estResultLATX ) )
print( summary( elas( estResultLATX ) ) )


cat( "\nLA: Elasticity formula of Eales + Unnevehr\n" )
ela <- aidsElas( estResultLA$coef, wMeans, method = "EU" )
print( ela )

cat( "\nLA: Elasticity formula of Green + Alston\n" )
ela <- aidsElas( estResultLA$coef, wMeans, pMeans, method = "GA" )
print( ela )

cat( "\nLA: Elasticity formula of Buse\n" )
ela <- aidsElas( estResultLA$coef, wMeans, pMeans, method = "B1" )
print( ela )

cat( "\nLA: Elasticity formula of Buse (alternative formula)\n" )
ela <- aidsElas( estResultLA$coef, wMeans, pMeans, method = "B2" )
print( ela )


############# Price indices ##############
options( digits = 5 )
cat( "\n************** Price indices **************\n" )
cat( "\nStone index\n" )
pxS <- aidsPx( "S", pNames, wNames, Blanciforti86 )
print( pxS )

cat( "\nStone index with lagged shares\n" )
pxSL <- aidsPx( "SL", pNames, wNames, Blanciforti86 )
print( pxSL )

cat( "\nPaasche index\n" )
pxP <- aidsPx( "P", pNames, wNames, Blanciforti86 )
print( pxP )

cat( "\nLaspeyres index\n" )
pxL <- aidsPx( "L", pNames, wNames, Blanciforti86 )
print( pxL )

cat( "\nTornqvist index\n" )
pxT <- aidsPx( "T", pNames, wNames, Blanciforti86 )
print( pxT )

cat( "\nTranslog index\n" )
pxTL <- aidsPx( "TL", pNames, wNames, Blanciforti86,
   coef = estResultLA$coef )
print( pxTL )

########### fitted values #################
options( digits = 3 )
fittedAIDS <- aidsCalc( pNames, "xFood", Blanciforti86[ -1, ],
   coef = estResultAIDS$coef )
print( fittedAIDS )
if( max( abs( fittedAIDS$shares[ !is.na( fittedAIDS$shares ) ] -
   estResultAIDS$wFitted ) ) > 1e-5 ) {
   stop( "fitted shares of AIDS are wrong" )
}
if( max( abs( fittedAIDS$quant[ !is.na( fittedAIDS$quant ) ] -
   estResultAIDS$qFitted ) ) > 1e-5 ) {
   stop( "fitted quantities of AIDS are wrong" )
}
fittedAIDSTX <- aidsCalc( pNames, "xFood", Blanciforti86[ -1, ],
   coef = estResultAIDSTX$coef )
print( fittedAIDSTX )
print( all.equal( fittedAIDS, fittedAIDSTX ) )

fittedLA <- aidsCalc( pNames, "xFood", Blanciforti86[ set, ],
   coef = estResultLA$coef, lnp = estResultLA$lnp )
print( fittedLA )
if( max( abs( fittedLA$shares[ -1, ] - estResultLA$wFitted[ setWo1, ] ),
      na.rm = TRUE ) > 1e-5 ) {
   stop( "fitted shares of LA-AIDS are wrong" )
}
if( max( abs( fittedLA$quant[ -1, ] - estResultLA$qFitted[ setWo1, ] ),
      na.rm = TRUE ) > 1e-5 ) {
   stop( "fitted quantities of LA-AIDS are wrong" )
}
fittedLATX <- aidsCalc( pNames, "xFood", Blanciforti86[ set, ],
   coef = estResultLATX$coef, lnp = estResultLATX$lnp )
print( fittedLATX )
print( all.equal( fittedLA, fittedLATX ) )

####### consistency ###################
consist <- aidsTestConsist( pNames, wNames, "xFood", Blanciforti86[ set, ],
   coef = estResultAIDS$coef )
print( consist )
class( consist ) <- NULL
print( consist )

