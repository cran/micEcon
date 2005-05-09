library( micEcon )
data( Mroz87 )
options( digits = 6 )

## Greene( 2003 ): example 22.8, page 786
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
greene <- heckit( wage ~ exper + I( exper^2 ) + educ + city,
   lfp ~ age + I( age^2 ) + faminc + kids + educ, Mroz87 )
print( summary( greene ) )
print( summary( greene$lm ) )
print( summary( greene$probit ) )
print( greene$sigma )
print( greene$rho )
print( greene$vcov )

## Wooldridge( 2003 ): example 17.5, page 590
data( Mroz87 )
wooldridge <- heckit( log( wage ) ~ educ + exper + I( exper^2 ),
   lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age + kids5 + kids618, Mroz87 )
print( summary( wooldridge ) )
print( summary( wooldridge$lm ) )
print( summary( wooldridge$probit ) )
print( wooldridge$sigma )
print( wooldridge$rho )
print( wooldridge$vcov )
