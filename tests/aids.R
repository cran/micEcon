library( micEcon )
data( Blanciforti86 )
options( digits = 6 )

pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
wMeans <- colMeans( Blanciforti86[ , wNames ] )
pMeans <- colMeans( Blanciforti86[ , pNames ] )

cat( paste( "\nRepeating the demand analysis of Blanciforti, Green",
   "& King (1986)\n" ) )
estResultLA <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86, method = "LA:SL", elaFormula = "Ch",
   maxiter = 1, rcovformula = 1, tol = 1e-7 )
print( estResultLA )

cat( paste( "\nRepeating the evaluation of different elasticity formulas",
   "of Green & Alston (1990): iterated AIDS\n" ) )
estResultAIDS <- aidsEst( pNames, wNames, "xFood",
   data = Blanciforti86[ -1, ], maxiter = 1,
   elaFormula = "AIDS", rcovformula=1, tol=1e-7,
   method = "MK:L" )
print( estResultAIDS )

########## Elasticities ###############
cat( "\nAIDS: Elasticities\n" )
wMeans <- colMeans( Blanciforti86[ , wNames ] )
ela <- aidsEla( estResultAIDS$coef, wMeans, pMeans, formula = "AIDS" )
print( ela )

cat( "\nLA: Elasticity formula of non-linear AIDS\n" )
wMeans <- colMeans( Blanciforti86[ , wNames ] )
ela <- aidsEla( estResultLA$coef, wMeans, pMeans, formula = "AIDS" )
print( ela )

cat( "\n********** Elasticities ***************" )
cat( "\nLA: Elasticity formula of Chalfant / Goddard\n" )
ela <- aidsEla( estResultLA$coef, wMeans, formula = "Ch" )
print( ela )

cat( "\nLA: Elasticity formula of Eales + Unnevehr\n" )
wMeans <- colMeans( Blanciforti86[ , wNames ] )
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
if( max( abs( fittedAIDS$shares - estResultAIDS$wFitted ) ) > 1e-5 ) {
   stop( "Fitted shares of AIDS are wrong." )
}
if( max( abs( fittedAIDS$quant - estResultAIDS$qFitted ) ) > 1e-5 ) {
   stop( "Fitted quantities of AIDS are wrong." )
}

fittedLA <- aidsCalc( pNames, "xFood", Blanciforti86,
   coef = estResultLA$coef, lnp = estResultLA$lnp )
print( fittedLA )
if( max( abs( fittedLA$shares[ -1, ] - estResultLA$wFitted[ -1, ] ) ) > 1e-5 ) {
   stop( "Fitted shares of LA-AIDS are wrong." )
}
if( max( abs( fittedLA$quant[ -1, ] - estResultLA$qFitted[ -1, ] ) ) > 1e-5 ) {
   stop( "Fitted quantities of LA-AIDS are wrong." )
}

####### consistency ###################
consist <- aidsTestConsist( pNames, wNames, "xFood", Blanciforti86,
   coef = estResultAIDS$coef )
print( consist )
class( consist ) <- NULL
print( consist )

