library( micEcon )
data( Mroz87 )
options( digits = 6 )

## Greene( 2003 ): example 22.8, page 786
Mroz87$kids  <- ( Mroz87$kids5 + Mroz87$kids618 > 0 )
greene <- heckit( lfp ~ age + I( age^2 ) + faminc + kids + educ,
   wage ~ exper + I( exper^2 ) + educ + city, Mroz87 )
print( greene )
print(summary( greene ))
print(summary( greene$lm ) )
print( summary( greene$probit ) )
print( greene$sigma )
print( greene$rho )
print( coef( greene ), digits = 5 )
print( coef( greene, part = "outcome" ), digits = 5 )
print( coef( summary( greene ) ), digits = 5 )
print( coef( summary( greene ), part = "outcome" ), digits = 5 )
print( vcov( greene ), digits = 5 )
print( vcov( greene, part = "outcome" ), digits = 5 )

## Wooldridge( 2003 ): example 17.5, page 590
data( Mroz87 )
wooldridge <- heckit( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, log( wage ) ~ educ + exper + I( exper^2 ), Mroz87 )
print( wooldridge )
print( summary( wooldridge ) )
print( summary( wooldridge$lm ) )
print( summary( wooldridge$probit ) )
print( wooldridge$sigma )
print( wooldridge$rho )
print( coef( wooldridge ), digits = 5 )
print( coef( wooldridge, part = "outcome" ), digits = 5 )
print( coef( summary( wooldridge ) ), digits = 5 )
print( coef( summary( wooldridge ), part = "outcome" ), digits = 5 )
print( vcov( wooldridge ), digits = 5 )
print( vcov( wooldridge, part = "outcome" ), digits = 5 )

## Tobit 5 Example from the selection paper
library(mvtnorm)
set.seed(0)
vc <- diag(3)
vc[lower.tri(vc)] <- c(0.9, 0.5, 0.1)
vc[upper.tri(vc)] <- vc[lower.tri(vc)]
eps <- rmvnorm(500, rep(0, 3), vc)
xs <- runif(500)
ys <- xs + eps[,1] > 0
xo1 <- runif(500)
yo1 <- xo1 + eps[,2]
xo2 <- runif(500)
yo2 <- xo2 + eps[,3]
heckit5test <- heckit( ys~xs, list( yo1 ~ xo1, yo2 ~ xo2 ) )
print( heckit5test )
print( summary( heckit5test ) )
print( coef( wooldridge ), digits = 5 )
print( coef( wooldridge, part = "outcome" ), digits = 5 )
print( coef( summary( wooldridge ) ), digits = 5 )
print( coef( summary( wooldridge ), part = "outcome" ), digits = 5 )
print( vcov( wooldridge ), digits = 5 )
print( vcov( wooldridge, part = "outcome" ), digits = 5 )
