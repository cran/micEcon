library( micEcon )

## Missong77
data( Missong77 )
## price indices for Missong77
cat( "\nLaspeyres Price Indices for Missong77\n" )
print( priceIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77 ) )

cat( "\nPaasche Price Indices for Missong77\n" )
print( priceIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77, "Paasche" ) )

cat( "\nFisher Price Indices for Missong77\n" )
print( priceIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77, "Fisher" ) )

## quantity indices for Missong77
cat( "\nLaspeyres Quantity Indices for Missong77\n" )
print( quantityIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77 ) )

cat( "\nPaasche Quantity Indices for Missong77\n" )
print( quantityIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77, "Paasche" ) )

cat( "\nFisher Quantity Indices for Missong77\n" )
print( quantityIndex( c( "p.beef", "p.veal", "p.pork" ),
   c( "q.beef", "q.veal", "q.pork" ), 1, Missong77, "Fisher" ) )


## Bleymueller251
data( Bleymueller251 )
## price indices for Bleymueller251
cat( "\nLaspeyres Price Indices for Bleymueller251\n" )
print( priceIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c( "q.A", "q.B", "q.C", "q.D" ),  1, Bleymueller251 ) )

cat( "\nPaasche Price Indices for Bleymueller251\n" )
print( priceIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c( "q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251, "Paasche" ) )

cat( "\nFisher Price Indices for Bleymueller251\n" )
print( priceIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251, "Fisher" ) )

## quantity indices for Bleymueller251
cat( "\nLaspeyres Quantity Indices for Bleymueller251\n" )
print( quantityIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251 ) )

cat( "\nPaasche Quantity Indices for Bleymueller251\n" )
print( quantityIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251, "Paasche" ) )

cat( "\nFisher Quantity Indices for Bleymueller251\n" )
print( quantityIndex( c( "p.A", "p.B", "p.C", "p.D" ),
   c("q.A", "q.B", "q.C", "q.D" ), 1, Bleymueller251, "Fisher" ) )


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
allObs <- rep( TRUE, nrow( Blanciforti86 ) )

## price indices for Blanciforti
cat( "\nLaspeyres Price Indices for Blanciforti\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86 ) )

cat( "\nPaasche Price Indices for Blanciforti\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Paasche" ) )

cat( "\nFisher Price Indices for Blanciforti\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Fisher" ) )

## quantity indices for Blanciforti
cat( "\nLaspeyres Quantity Indices for Blanciforti\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86 ) )

cat( "\nPaasche Quantity Indices for Blanciforti\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Paasche" ) )

cat( "\nFisher Quantity Indices for Blanciforti\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Fisher" ) )

## price indices for Blanciforti and base=mean
cat( "\nLaspeyres Price Indices for Blanciforti and base=mean\n" )
print( priceIndex( pNames, qNames, allObs, Blanciforti86 ) )

cat( "\nPaasche Price Indices for Blanciforti and base=mean\n" )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Paasche" ) )

cat( "\nFisher Price Indices for Blanciforti and base=mean\n" )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Fisher" ) )

## quantity indices for Blanciforti and base=mean
cat( "\nLaspeyres Quantity Indices for Blanciforti and base=mean\n" )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86 ) )

cat( "\nPaasche Quantity Indices for Blanciforti and base=mean\n" )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Paasche" ) )

cat( "\nFisher Quantity Indices for Blanciforti and base=mean\n" )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Fisher" ) )


## Blanciforti with some NA prices
## manipulating data of Blanciforti
for( i in 1:4 ) {
   Blanciforti86[[ pNames[ i ] ]][ c( 2, i * 4, i * 8 ) ] <- NA
}

## price indices for Blanciforti with some NA prices
cat( "\nLaspeyres Price Indices for Blanciforti with some NA prices\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86 ) )

cat( "\nPaasche Price Indices for Blanciforti with some NA prices\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Paasche" ) )

cat( "\nFisher Price Indices for Blanciforti with some NA prices\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Fisher" ) )

## quantity indices for Blanciforti with some NA prices
cat( "\nLaspeyres Quantity Indices for Blanciforti with some NA prices\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86 ) )

cat( "\nPaasche Quantity Indices for Blanciforti with some NA prices\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Paasche" ) )

cat( "\nFisher Quantity Indices for Blanciforti with some NA prices\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Fisher" ) )

## price indices for Blanciforti with some NA prices and na.rm=TRUE
cat( "\nLaspeyres Price Indices for Blanciforti with some NA prices and na.rm=TRUE\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, na.rm = TRUE ) )

cat( "\nPaasche Price Indices for Blanciforti with some NA prices and na.rm=TRUE\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Paasche", na.rm = TRUE ) )

cat( "\nFisher Price Indices for Blanciforti with some NA prices and na.rm=TRUE\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Fisher", na.rm = TRUE ) )

## quantity indices for Blanciforti with some NA prices and na.rm=TRUE
cat( "\nLaspeyres Quantity Indices for Blanciforti with some NA prices and na.rm=TRUE\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, na.rm = TRUE ) )

cat( "\nPaasche Quantity Indices for Blanciforti with some NA prices and na.rm=TRUE\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Paasche", na.rm = TRUE ) )

cat( "\nFisher Quantity Indices for Blanciforti with some NA prices and na.rm=TRUE\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Fisher", na.rm = TRUE ) )

## price indices for Blanciforti with some NA prices and base=mean
cat( "\nLaspeyres Price Indices for Blanciforti with some NA prices and base=mean\n" )
print( priceIndex( pNames, qNames, 16, Blanciforti86 ) )

cat( "\nPaasche Price Indices for Blanciforti with some NA prices and base=mean\n" )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Paasche" ) )

cat( "\nFisher Price Indices for Blanciforti with some NA prices and base=mean\n" )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Fisher" ) )

## quantity indices for Blanciforti with some NA prices and base=mean
cat( "\nLaspeyres Quantity Indices for Blanciforti with some NA prices and base=mean\n" )
print( quantityIndex( pNames, qNames, 16, Blanciforti86 ) )

cat( "\nPaasche Quantity Indices for Blanciforti with some NA prices and base=mean\n" )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Paasche" ) )

cat( "\nFisher Quantity Indices for Blanciforti with some NA prices and base=mean\n" )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Fisher" ) )

## price indices for Blanciforti with some NA prices and na.rm=TRUE and base=mean
cat( paste( "\nLaspeyres Price Indices for Blanciforti with some NA prices",
   "and na.rm=TRUE and base=mean\n" ) )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, na.rm = TRUE ) )

cat( paste( "\nPaasche Price Indices for Blanciforti with some NA prices",
   "and na.rm=TRUE and base=mean\n" ) )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Paasche", na.rm = TRUE ) )

cat( paste( "\nFisher Price Indices for Blanciforti with some NA prices",
   "and na.rm=TRUE and base=mean\n" ) )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Fisher", na.rm = TRUE ) )

## quantity indices for Blanciforti with some NA prices and na.rm=TRUE and base=mean
cat( paste( "\nLaspeyres Quantity Indices for Blanciforti with some NA prices",
   "and na.rm=TRUE and base=mean\n" ) )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, na.rm = TRUE ) )

cat( paste( "\nPaasche Quantity Indices for Blanciforti with some NA prices",
   "and na.rm=TRUE and base=mean\n" ) )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Paasche", na.rm = TRUE ) )

cat( paste( "\nFisher Quantity Indices for Blanciforti with some NA prices",
   "and na.rm=TRUE and base=mean\n" ) )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Fisher", na.rm = TRUE ) )


## Blanciforti with some NA prices and quantities
## manipulating data of Blanciforti
for( i in 1:4 ) {
   Blanciforti86[[ qNames[ i ] ]][ c( 2, ( i + 1 ) * 4, i * 8 ) ] <- NA
}

## price indices for Blanciforti with some NA prices and quantities
cat( "\nLaspeyres Price Indices for Blanciforti with NA prices and quantities\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86 ) )

cat( "\nPaasche Price Indices for Blanciforti with NA prices and quantities\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Paasche" ) )

cat( "\nFisher Price Indices for Blanciforti with NA prices and quantities\n" )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Fisher" ) )

## quantity indices for Blanciforti with some NA prices and quantities
cat( "\nLaspeyres Quantity Indices for Blanciforti with NA prices and quantities\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86 ) )

cat( "\nPaasche Quantity Indices for Blanciforti with NA prices and quantities\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Paasche" ) )

cat( "\nFisher Quantity Indices for Blanciforti with NA prices and quantities\n" )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Fisher" ) )

## price indices for Blanciforti with some NA prices and quantities and na.rm=TRUE
cat( paste( "\nLaspeyres Price Indices for Blanciforti with NA prices and",
   "quantities and na.rm=TRUE\n" ) )
print( priceIndex( pNames, qNames, 1, Blanciforti86, na.rm = TRUE ) )

cat( paste( "\nPaasche Price Indices for Blanciforti with NA prices and",
   "quantities and na.rm=TRUE\n" ) )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Paasche", na.rm = TRUE ) )

cat( paste( "\nFisher Price Indices for Blanciforti with NA prices and",
   "quantities and na.rm=TRUE\n" ) )
print( priceIndex( pNames, qNames, 1, Blanciforti86, "Fisher", na.rm = TRUE ) )

## quantity indices for Blanciforti with some NA prices and quantities and na.rm=TRUE
cat( paste( "\nLaspeyres Quantity Indices for Blanciforti with NA prices and",
   "quantities and na.rm=TRUE\n" ) )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, na.rm = TRUE ) )

cat( paste( "\nPaasche Quantity Indices for Blanciforti with NA prices and",
   "quantities and na.rm=TRUE\n" ) )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Paasche", na.rm = TRUE ) )

cat( paste( "\nFisher Quantity Indices for Blanciforti with some NA prices and",
   "na.rm=TRUE\n" ) )
print( quantityIndex( pNames, qNames, 1, Blanciforti86, "Fisher", na.rm = TRUE ) )

## price indices for Blanciforti with some NA prices and quantities and base=mean
cat( paste( "\nLaspeyres Price Indices for Blanciforti with NA prices and",
   "quantities and base=mean\n" ) )
print( priceIndex( pNames, qNames, 16, Blanciforti86 ) )

cat( paste( "\nPaasche Price Indices for Blanciforti with NA prices and",
   "quantities and base=mean\n" ) )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Paasche" ) )

cat( paste( "\nFisher Price Indices for Blanciforti with NA prices and",
   "quantities and base=mean\n" ) )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Fisher" ) )

## quantity indices for Blanciforti with some NA prices and quantities and base=mean
cat( paste( "\nLaspeyres Quantity Indices for Blanciforti with NA prices and",
   "quantities and base=mean\n" ) )
print( quantityIndex( pNames, qNames, 16, Blanciforti86 ) )

cat( paste( "\nPaasche Quantity Indices for Blanciforti with NA prices and",
   "quantities and base=mean\n" ) )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Paasche" ) )

cat( paste( "\nFisher Quantity Indices for Blanciforti with NA prices and",
   "quantities and base=mean\n" ) )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Fisher" ) )

## price indices for Blanciforti with some NA prices and quantities and na.rm=TRUE and base=mean
cat( paste( "\nLaspeyres Price Indices for Blanciforti with NA prices and quantities",
   "and na.rm=TRUE and base=mean\n" ) )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, na.rm = TRUE ) )

cat( paste( "\nPaasche Price Indices for Blanciforti with NA prices and quantities",
   "and na.rm=TRUE and base=mean\n" ) )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Paasche", na.rm = TRUE ) )

cat( paste( "\nFisher Price Indices for Blanciforti with NA prices and quantities",
   "and na.rm=TRUE and base=mean\n" ) )
print( priceIndex( pNames, qNames, allObs, Blanciforti86, "Fisher", na.rm = TRUE ) )

## quantity indices for Blanciforti with some NA prices and quantities and na.rm=TRUE and base=mean
cat( paste( "\nLaspeyres Quantity Indices for Blanciforti with NA prices and quantities",
   "and na.rm=TRUE and base=mean\n" ) )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, na.rm = TRUE ) )

cat( paste( "\nPaasche Quantity Indices for Blanciforti with NA prices and quantities",
   "and na.rm=TRUE and base=mean\n" ) )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Paasche", na.rm = TRUE ) )

cat( paste( "\nFisher Quantity Indices for Blanciforti with NA prices and quantities",
   "and na.rm=TRUE and base=mean\n" ) )
print( quantityIndex( pNames, qNames, allObs, Blanciforti86, "Fisher", na.rm = TRUE ) )

