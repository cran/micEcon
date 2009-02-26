library( micEcon )


## matrix
m <- matrix( 1:24, nrow = 6, ncol = 4 )

cm1 <- colMedians( m )
print( cm1 )

rm1 <- rowMedians( m )
print( rm1 )

all.equal( cm1, rowMedians( t( m ) ) )
all.equal( rm1, colMedians( t( m ) ) )


## data.frame
data( germanFarms )

cm2 <- colMedians( germanFarms[ , -1 ] )
print( cm2 )

rm2 <- rowMedians( germanFarms[ , -1 ] )
print( rm2 )

all.equal( cm2, rowMedians( t( germanFarms[ , -1 ] ) ) )
all.equal( rm2, colMedians( t( germanFarms[ , -1 ] ) ) )
