library( micEcon )

## function for printIndexing indices
printIndices <- function( what, ... ) {
   for( i in c( "Laspeyres", "Paasche", "Fisher" ) ) {
      cat( "\n", i, "\n" )
      if( what == "p" ) {
         index <- priceIndex( ..., method = i, weights = TRUE )
      } else {
         index <- quantityIndex( ..., method = i, weights = TRUE )
      }
      print( index )
      testRowSums <- rowSums( attributes( index )$weights[ !is.na( index ), ] )
      names( testRowSums ) <- NULL
      if( all.equal( testRowSums,
         rep( 1, sum( !is.na( index ) ) ) ) != TRUE ) {
         cat( "\nrowSums are not equal to one!!!\n\n" )
      }
   }
}

## Missong03E7.7
data( Missong03E7.7 )
## price indices for Missong03E7.7
printIndices( "p",  c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong03E7.7 )

## quantity indices for Missong03E7.7
printIndices( "q",  c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong03E7.7 )


## Bleymueller79E25.1
data( Bleymueller79E25.1 )
## price indices for Bleymueller79E25.1
printIndices( "p",  c( "p.A", "p.B", "p.C", "p.D" ),
   c( "q.A", "q.B", "q.C", "q.D" ),  1, Bleymueller79E25.1 )

## quantity indices for Bleymueller79E25.1
printIndices( "q",  c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller79E25.1 )


## Blanciforti
data( Blanciforti86 )
## preparing data of Blanciforti
pNames <- c( "pFood1", "pFood2", "pFood3", "pFood4" )
qNames <- c( "qFood1", "qFood2", "qFood3", "qFood4" )
wNames <- c( "wFood1", "wFood2", "wFood3", "wFood4" )
for( i in 1:4 ) {
   Blanciforti86[[ qNames[ i ] ]] <- Blanciforti86[[ wNames[ i ] ]] *
      Blanciforti86[[ "xFood" ]] / Blanciforti86[[ pNames[ i ] ]]
}
allObs <- !is.na( Blanciforti86$pFood1 )

## price indices for Blanciforti
printIndices( "p",  pNames, qNames, 1, Blanciforti86 )

## quantity indices for Blanciforti
printIndices( "q",  pNames, qNames, 1, Blanciforti86 )

## price indices for Blanciforti and base=mean
printIndices( "p",  pNames, qNames, allObs, Blanciforti86 )

## quantity indices for Blanciforti and base=mean
printIndices( "q",  pNames, qNames, allObs, Blanciforti86 )


## Blanciforti with some NA prices
## manipulating data of Blanciforti
for( i in 1:4 ) {
   Blanciforti86[[ pNames[ i ] ]][ c( 2, i * 4, i * 8 ) ] <- NA
}

## price indices for Blanciforti with some NA prices
printIndices( "p",  pNames, qNames, 1, Blanciforti86 )

## quantity indices for Blanciforti with some NA prices
printIndices( "q",  pNames, qNames, 1, Blanciforti86 )

## price indices for Blanciforti with some NA prices and na.rm=TRUE
printIndices( "p",  pNames, qNames, 1, Blanciforti86, na.rm = TRUE )

## quantity indices for Blanciforti with some NA prices and na.rm=TRUE
printIndices( "q",  pNames, qNames, 1, Blanciforti86, na.rm = TRUE )

## price indices for Blanciforti with some NA prices and base=mean
printIndices( "p",  pNames, qNames, 16, Blanciforti86 )

## quantity indices for Blanciforti with some NA prices and base=mean
printIndices( "q",  pNames, qNames, allObs, Blanciforti86 )

## price indices for Blanciforti with some NA prices and na.rm=TRUE and base=mean
printIndices( "p",  pNames, qNames, allObs, Blanciforti86, na.rm = TRUE )

## quantity indices for Blanciforti with some NA prices and na.rm=TRUE and base=mean
printIndices( "q",  pNames, qNames, allObs, Blanciforti86, na.rm = TRUE )


## Blanciforti with some NA prices and quantities
## manipulating data of Blanciforti
for( i in 1:4 ) {
   Blanciforti86[[ qNames[ i ] ]][ c( 2, ( i + 1 ) * 4, i * 8 ) ] <- NA
}

## price indices for Blanciforti with some NA prices and quantities
printIndices( "p",  pNames, qNames, 1, Blanciforti86 )

## quantity indices for Blanciforti with some NA prices and quantities
printIndices( "q",  pNames, qNames, 1, Blanciforti86 )

## price indices for Blanciforti with some NA prices and quantities and na.rm=TRUE
printIndices( "p",  pNames, qNames, 1, Blanciforti86, na.rm = TRUE )

## quantity indices for Blanciforti with some NA prices and quantities and na.rm=TRUE
printIndices( "q",  pNames, qNames, 1, Blanciforti86, na.rm = TRUE )

## price indices for Blanciforti with some NA prices and quantities and base=mean
printIndices( "p",  pNames, qNames, 16, Blanciforti86 )

## quantity indices for Blanciforti with some NA prices and quantities and base=mean
printIndices( "q",  pNames, qNames, 16, Blanciforti86 )

## price indices for Blanciforti with some NA prices and quantities and na.rm=TRUE and base=mean
printIndices( "p",  pNames, qNames, allObs, Blanciforti86, na.rm = TRUE )

## quantity indices for Blanciforti with some NA prices and quantities and na.rm=TRUE and base=mean
printIndices( "q",  pNames, qNames, allObs, Blanciforti86, na.rm = TRUE )
