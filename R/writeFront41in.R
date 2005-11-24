writeFront41in <- function( data, crossSectionName, timePeriodName,
   yName, xNames = NULL, zNames = NULL,
   translog = FALSE, quadHalf = TRUE,
   functionType = 1, modelType = 1, logDepVar = TRUE, mu = FALSE, eta = FALSE,
   insFile = "front41.ins", dtaFile = sub( "\.ins$", ".dta", insFile ),
   outFile = sub( "\.ins$", ".out", insFile ) ) {

   checkNames( c( crossSectionName, timePeriodName, yName, xNames, zNames ),
      names( data ) )

   if( !modelType %in% c( 1, 2 ) ) {
      stop( "argument 'modelType' must be either 1 or 2" )
   }
   if( !functionType %in% c( 1, 2 ) ) {
      stop( "argument 'functionType' must be either 1 or 2" )
   }
   if( !is.logical( logDepVar ) ) {
      stop( "argument 'logDepVar' must be logical" )
   }
   if( !is.logical( mu ) ) {
      stop( "argument 'mu' must be logical" )
   }
   if( modelType == 1 ) {
      if( !is.logical( eta ) ) {
         stop( "argument 'eta' must be logical" )
      }
   }


   nCrossSection <- max( data[[ crossSectionName ]] )
   nTimePeriods  <- max( data[[ timePeriodName ]] )
   nTotalObs     <- nrow( data )
   nXvars        <- length( xNames )
   nXtotal       <- ifelse( translog, nXvars + nXvars * ( nXvars + 1 ) / 2,
                            nXvars )
   nZvars        <- length( zNames )

   if( modelType == 2 ) {
      eta <- nZvars
   } else {
      eta <- ifelse( eta, "y", "n" )
   }

   commentRow <- max( 16, nchar( dtaFile ) + 1 )

   cat( modelType, rep( " ", commentRow - 1 ),
      "1=ERROR COMPONENTS MODEL, 2=TE EFFECTS MODEL\n",
      file = insFile, sep = "" )
   cat( dtaFile, rep( " ", commentRow - nchar( dtaFile ) ),
      "DATA FILE NAME\n", file = insFile, append = TRUE, sep = "" )
   cat( outFile, rep( " ", commentRow - nchar( outFile ) ),
      "OUTPUT FILE NAME\n", file = insFile, append = TRUE, sep = "" )
   cat( functionType, rep( " ", commentRow - 1 ),
      "1=PRODUCTION FUNCTION, 2=COST FUNCTION\n",
      file = insFile, append = TRUE, sep = "" )
   cat( ifelse( logDepVar, "y", "n" ), rep( " ", commentRow - 1 ),
      "LOGGED DEPENDENT VARIABLE (Y/N)\n",
      file = insFile, append = TRUE, sep = "" )
   cat( nCrossSection,
      rep( " ", commentRow - nchar( as.character( nCrossSection ) ) ),
      "NUMBER OF CROSS-SECTIONS\n",
      file = insFile, append = TRUE, sep = "" )
   cat( nTimePeriods,
      rep( " ", commentRow - nchar( as.character( nTimePeriods ) ) ),
      "NUMBER OF TIME PERIODS\n",
      file = insFile, append = TRUE, sep = "" )
   cat( nTotalObs,
      rep( " ", commentRow - nchar( as.character( nTotalObs ) ) ),
      "NUMBER OF OBSERVATIONS IN TOTAL\n",
      file = insFile, append = TRUE, sep = "" )
   cat( nXtotal,
      rep( " ", commentRow - nchar( as.character( nXtotal ) ) ),
      "NUMBER OF REGRESSOR VARIABLES (Xs)\n",
      file = insFile, append = TRUE, sep = "" )
   cat( ifelse( mu, "y", "n" ), rep( " ", commentRow - 1 ),
      "MU (Y/N) [OR DELTA0 (Y/N) IF USING TE EFFECTS MODEL]\n",
      file = insFile, append = TRUE, sep = "" )
   cat( eta, rep( " ", commentRow - nchar( as.character( eta ) ) ),
      "ETA (Y/N) [OR NUMBER OF TE EFFECTS REGRESSORS (Zs)]\n",
      file = insFile, append = TRUE, sep = "" )
   cat( "n", rep( " ", commentRow - 1 ),
      "STARTING VALUES (Y/N)\n",
      file = insFile, append = TRUE, sep = "" )

   dataTable <- cbind( data[[ crossSectionName ]], data[[ timePeriodName ]],
      data[[ yName ]] )

   if( nXvars > 0 ) {
      for( i in 1:nXvars ) {
         dataTable <- cbind( dataTable, data[[ xNames[ i ] ]] )
      }
      if( translog ) {
         for( i in 1:nXvars ) {
            for( j in i:nXvars ) {
               dataTable <- cbind( dataTable,
                  ifelse( i == j, 1 , 2 ) * ifelse( quadHalf, 0.5, 1 ) *
                  data[[ xNames[ i ] ]] * data[[ xNames[ j ] ]] )
            }
         }
      }
   }
   if( nZvars > 0 ) {
      for( i in 1:nZvars ) {
         dataTable <- cbind( dataTable, data[[ zNames[ i ] ]] )
      }
   }
   write.table( dataTable, file = dtaFile, row.names = FALSE,
      col.names = FALSE, sep = "\t" )
}