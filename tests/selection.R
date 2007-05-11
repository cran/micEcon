library(micEcon)
library(mvtnorm)
options(digits=6)
N <- 1500
NNA <- 5
vc <- diag(3)
vc[lower.tri(vc)] <- c(0.9, 0.5, 0.6)
vc[upper.tri(vc)] <- vc[lower.tri(vc)]
set.seed(1)
## ------- Tobit-5 example ---------
eps <- rmvnorm(N, rep(0, 3), vc)
xs <- runif(N)
ys <- xs + eps[,1] > 0
xo1 <- runif(N)
yo1 <- xo1 + eps[,2]
xo2 <- runif(N)
yo2 <- xo2 + eps[,3]
## Put some NA-s into the data
ys[sample(N, NNA)] <- NA
xs[sample(N, NNA)] <- NA
xo1[sample(N, NNA)] <- NA
xo2[sample(N, NNA)] <- NA
yo1[sample(N, NNA)] <- NA
yo2[sample(N, NNA)] <- NA
testTobit5TwoStep <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), method="2step")
print( testTobit5TwoStep )
print( summary( testTobit5TwoStep ) )
print( coef( testTobit5TwoStep ), digits = 5 )
print( coef( testTobit5TwoStep, part = "outcome" ), digits = 5 )
print( coef( summary( testTobit5TwoStep ) ), digits = 5 )
print( coef( summary( testTobit5TwoStep ), part = "outcome" ), digits = 5 )
print( vcov( testTobit5TwoStep ), digits = 5 )
print( vcov( testTobit5TwoStep, part = "outcome" ), digits = 5 )

testTobit5Ml <- selection(ys~xs, list(yo1 ~ xo1, yo2 ~ xo2), method="ml")
print( testTobit5Ml )
print( summary( testTobit5Ml ) )
print( coef( testTobit5Ml ), digits = 5 )
print( coef( testTobit5Ml, part = "outcome" ), digits = 5 )
print( coef( summary( testTobit5Ml ) ), digits = 5 )
print( coef( summary( testTobit5Ml ), part = "outcome" ), digits = 5 )
print( vcov( testTobit5Ml ), digits = 5 )
print( vcov( testTobit5Ml, part = "outcome" ), digits = 5 )

## ------- Tobit-2 exmple -----------
vc <- diag(2)
vc[2,1] <- vc[1,2] <- -0.7
eps <- rmvnorm(N, rep(0, 2), vc)
xs <- runif(N)
ys <- xs + eps[,1] > 0
xo <- runif(N)
yo <- (xo + eps[,2])*(ys > 0)
xs[sample(N, NNA)] <- NA
ys[sample(N, NNA)] <- NA
xo[sample(N, NNA)] <- NA
yo[sample(N, NNA)] <- NA
testTobit2TwoStep <- selection(ys~xs, yo ~xo, method="2step")
print( testTobit2TwoStep )
print( summary( testTobit2TwoStep ) )
print( coef( testTobit2TwoStep ), digits = 5 )
print( coef( testTobit2TwoStep, part = "outcome" ), digits = 5 )
print( coef( summary( testTobit2TwoStep ) ), digits = 5 )
print( coef( summary( testTobit2TwoStep ), part = "outcome" ), digits = 5 )
print( vcov( testTobit2TwoStep ), digits = 5 )
print( vcov( testTobit2TwoStep, part = "outcome" ), digits = 5 )

testTobit2Ml <- selection(ys~xs, yo ~xo, method="ml")
print( testTobit2Ml )
print( summary( testTobit2Ml ) )
print( coef( testTobit2Ml ), digits = 5 )
print( coef( testTobit2Ml, part = "outcome" ), digits = 5 )
print( coef( summary( testTobit2Ml ) ), digits = 5 )
print( coef( summary( testTobit2Ml ), part = "outcome" ), digits = 5 )
print( vcov( testTobit2Ml ), digits = 5 )
print( vcov( testTobit2Ml, part = "outcome" ), digits = 5 )
