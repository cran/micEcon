library( micEcon )
data( Mroz87 )
options( digits = 6 )

## Greene( 2003 ): example 22.8, page 786
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
greene <- heckit( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
print( summary( greene ) )
print( summary( greene$lm ) )
print( summary( greene$probit ) )
print( greene$sigma )
print( greene$rho )
print( greene$vcov, digits = 5 )

## Wooldridge( 2003 ): example 17.5, page 590
data( Mroz87 )
wooldridge <- heckit( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87 )
print( summary( wooldridge ) )
print( summary( wooldridge$lm ) )
print( summary( wooldridge$probit ) )
print( wooldridge$sigma )
print( wooldridge$rho )
print( wooldridge$vcov, digits = 5 )
