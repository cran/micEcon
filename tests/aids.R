library( micEcon )
data( Blanciforti86 )
options( digits = 6 )

set <- !is.na( Blanciforti86$pFood1 )
setWo1 <- set & rownames( Blanciforti86 ) != 1947
pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
pMeans <- colMeans( Blanciforti86[ set, pNames ] )

cat( paste( "\nRepeating the demand analysis of Blanciforti, Green",
   "& King (1986)\n" ) )
estResultLA <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "LA:SL", elaFormula = "Ch",
   maxiter = 1, rcovformula = 1, tol = 1e-7 )
print( estResultLA )
# imposing restrictions via TX
estResultLATX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ set, ], method = "LA:SL", elaFormula = "Ch",
   maxiter = 1, rcovformula = 1, tol = 1e-7, TX = TRUE )
print( estResultLATX )
estResultLATX$est$bt <- NULL
estResultLATX$est$btcov <- NULL
estResultLATX$est$x <- NULL
estResultLATX$est$TX <- NULL
estResultLATX$est$q.restr <- NULL
estResultLA$est$x <- NULL
estResultLA$est$R.restr <- NULL
estResultLA$est$q.restr <- NULL
print( all.equal( estResultLA, estResultLATX ) )

## only homogeneity (no symmetry imposed)
estResultLAhom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL", elaFormula = "Ch",
   maxiter = 1, rcovformula = 1, tol = 1e-7 )
print( estResultLAhom )
# imposing restrictions via TX
estResultLAhomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL", elaFormula = "Ch",
   maxiter = 1, rcovformula = 1, tol = 1e-7, TX = TRUE )
print( estResultLAhomTX )
estResultLAhomTX$est$bt <- NULL
estResultLAhomTX$est$btcov <- NULL
estResultLAhomTX$est$x <- NULL
estResultLAhomTX$est$TX <- NULL
estResultLAhomTX$est$q.restr <- NULL
estResultLAhom$est$x <- NULL
estResultLAhom$est$R.restr <- NULL
estResultLAhom$est$q.restr <- NULL
print( all.equal( estResultLAhom, estResultLAhomTX ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultLAunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL", elaFormula = "Ch",
   maxiter = 1, rcovformula = 1, tol = 1e-7 )
print( estResultLAunr )
# imposing restrictions via TX
estResultLAunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ set, ], method = "LA:SL", elaFormula = "Ch",
   maxiter = 1, rcovformula = 1, tol = 1e-7, TX = TRUE )
print( estResultLAunrTX )
estResultLAunrTX$est$bt <- NULL
estResultLAunrTX$est$btcov <- NULL
estResultLAunrTX$est$x <- NULL
estResultLAunrTX$est$TX <- NULL
estResultLAunrTX$est$q.restr <- NULL
estResultLAunr$est$x <- NULL
estResultLAunr$est$R.restr <- NULL
estResultLAunr$est$q.restr <- NULL
print( all.equal( estResultLAunr, estResultLAunrTX ) )

#####################################################
cat( paste( "\nRepeating the evaluation of different elasticity formulas",
   "of Green & Alston (1990): iterated AIDS\n" ) )
estResultAIDS <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   elaFormula = "AIDS", rcovformula=1, tol=1e-7,
   method = "IL:L" )
print( estResultAIDS )
# imposing restrictions via TX
estResultAIDSTX <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   elaFormula = "AIDS", rcovformula=1, tol=1e-7,
   method = "IL:L", TX = TRUE )
print( estResultAIDSTX )
estResultAIDSTX$est$bt <- NULL
estResultAIDSTX$est$btcov <- NULL
estResultAIDSTX$est$x <- NULL
estResultAIDSTX$est$TX <- NULL
estResultAIDSTX$est$q.restr <- NULL
estResultAIDS$est$x <- NULL
estResultAIDS$est$R.restr <- NULL
estResultAIDS$est$q.restr <- NULL
print( all.equal( estResultAIDS, estResultAIDSTX ) )

## only homogeneity (no symmetry imposed)
estResultAIDShom <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   elaFormula = "AIDS", rcovformula=1, tol=1e-7,
   method = "IL:L" )
print( estResultAIDShom )
# imposing restrictions via TX
estResultAIDShomTX <- aidsEst( pNames, wNames, "xFood", sym = FALSE,
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   elaFormula = "AIDS", rcovformula=1, tol=1e-7,
   method = "IL:L", TX = TRUE )
print( estResultAIDShomTX )
estResultAIDShomTX$est$bt <- NULL
estResultAIDShomTX$est$btcov <- NULL
estResultAIDShomTX$est$x <- NULL
estResultAIDShomTX$est$TX <- NULL
estResultAIDShomTX$est$q.restr <- NULL
estResultAIDShom$est$x <- NULL
estResultAIDShom$est$R.restr <- NULL
estResultAIDShom$est$q.restr <- NULL
print( all.equal( estResultAIDShom, estResultAIDShomTX ) )

## unrestricted (no homogeneity and no symmetry imposed)
estResultAIDSunr <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   elaFormula = "AIDS", rcovformula=1, tol=1e-7,
   method = "IL:L" )
print( estResultAIDSunr )
# imposing restrictions via TX
estResultAIDSunrTX <- aidsEst( pNames, wNames, "xFood", hom = FALSE, sym = FALSE,
   data = Blanciforti86[ setWo1, ], maxiter = 1,
   elaFormula = "AIDS", rcovformula=1, tol=1e-7,
   method = "IL:L", TX = TRUE )
print( estResultAIDSunrTX )
estResultAIDSunrTX$est$bt <- NULL
estResultAIDSunrTX$est$btcov <- NULL
estResultAIDSunrTX$est$x <- NULL
estResultAIDSunrTX$est$TX <- NULL
estResultAIDSunrTX$est$q.restr <- NULL
estResultAIDSunr$est$x <- NULL
estResultAIDSunr$est$R.restr <- NULL
estResultAIDSunr$est$q.restr <- NULL
print( all.equal( estResultAIDSunr, estResultAIDSunrTX ) )

## with NAs
estResultLaSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, maxiter = 1, elaFormula = "AIDS",
   rcovformula=1, tol=1e-7, method = "LA:S" )
print( estResultLaSNa )

estResultLaSlNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, maxiter = 1, elaFormula = "AIDS",
   rcovformula=1, tol=1e-7, method = "LA:SL" )
print( estResultLaSlNa )

estResultLaLNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, maxiter = 1, elaFormula = "AIDS",
   rcovformula=1, tol=1e-7, method = "LA:L" )
print( estResultLaLNa )

estResultAIDSNa <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, maxiter = 1, elaFormula = "AIDS",
   rcovformula=1, tol=1e-7, method = "IL:L" )
print( estResultAIDSNa )


########## Elasticities ###############
cat( "\nAIDS: Elasticities\n" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
ela <- aidsEla( estResultAIDS$coef, wMeans, pMeans, formula = "AIDS" )
print( ela )
elaTX <- aidsEla( estResultAIDSTX$coef, wMeans, pMeans, formula = "AIDS" )
print( elaTX )
print( all.equal( ela, elaTX ) )


cat( "\nLA: Elasticity formula of non-linear AIDS\n" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
ela <- aidsEla( estResultLA$coef, wMeans, pMeans, formula = "AIDS" )
print( ela )
elaTX <- aidsEla( estResultLATX$coef, wMeans, pMeans, formula = "AIDS" )
print( elaTX )
print( all.equal( ela, elaTX ) )

cat( "\n********** Elasticities ***************" )
cat( "\nLA: Elasticity formula of Chalfant / Goddard\n" )
ela <- aidsEla( estResultLA$coef, wMeans, formula = "Ch" )
print( ela )

cat( "\nLA: Elasticity formula of Eales + Unnevehr\n" )
wMeans <- colMeans( Blanciforti86[ set, wNames ] )
ela <- aidsEla( estResultLA$coef, wMeans, formula = "EU" )
print( ela )


############# Price indices ##############
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

