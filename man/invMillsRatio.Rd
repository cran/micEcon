\name{invMillsRatio}
\alias{invMillsRatio}
\title{Inverse Mill's Ratio of probit models}

\description{
   Calculates the 'Inverse Mill's Ratios' of univariate and bivariate
   probit models.
}

\usage{
   invMillsRatio( x, all = FALSE )
}

\arguments{
   \item{x}{probit model estimated by \code{\link{glm}} or
      \code{\link[VGAM]{vglm}}.}
   \item{all}{a logical value indicating whether the inverse Mill's Ratios
      should be calculated for all observations.}
}

\details{
   The formula to calculate the inverse Mill's ratios for univariate probit
   models is taken from Greene (2003, p. 785),
   whereas the formulas for bivariate probit models are derived in
   Henning and Henningsen (2005).
}

\value{
   A data frame that contains the Inverse Mill's Ratios (IMR) and the delta
   values (see Greene, 2003, p. 784).

   If a univariate probit estimation is provided, the variables
   \code{IMR1} and \code{IMR0} are the Inverse Mill's Ratios to correct
   for a sample selection bias of y = 1 and y = 0, respectively.
   Accordingly, 'delta1' and 'delta0' are the corresponding delta values.

   If a bivariate probit estimation is provided, the variables
   \code{IMRa1}, \code{IMRa0}, \code{IMRb1}, and \code{IMRb0} are the
   Inverse Mills Ratios to correct for a sample selection bias
   of y = 1 and y = 0 in equations 'a' and 'b', respectively.
   Accordingly, 'deltaa1', 'deltaa0', 'deltab1' and 'deltab0' are the
   corresponding delta values.
}

\references{
  Greene, W. H. (2003)
  \emph{Econometric Analysis, Fifth Edition}, Prentice Hall.

  Henning, C.H.C.A and A. Henningsen (2005)
  Modeling Price Response of Farm Households in Imperfect
  Labor Markets in Poland:
  Incorporating Transaction Costs and Heterogeneity into a Farm
  Household Approach.
  Unpublished, University of Kiel, Germany.
}

\author{Arne Henningsen \email{ahenningsen@agric-econ.uni-kiel.de}}

\examples{
## Wooldridge( 2003 ): example 17.5, page 590
data(Mroz87)
myProbit <- glm( lfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, family = binomial( link = "probit" ), data=Mroz87 )
Mroz87$IMR <- invMillsRatio( myProbit )$IMR1
myHeckit <- lm( log( wage ) ~ educ + exper + I( exper^2 ) + IMR,
   data = Mroz87[ Mroz87$lfp == 1, ] )

# using NO labor force participation as endogenous variable
Mroz87$nolfp <- 1 - Mroz87$lfp
myProbit2 <- glm( nolfp ~ nwifeinc + educ + exper + I( exper^2 ) + age +
   kids5 + kids618, family = binomial( link = "probit" ), data=Mroz87 )
all.equal( invMillsRatio( myProbit )$IMR1, invMillsRatio( myProbit2 )$IMR0 )
   # should be true

# example for bivariate probit
\dontrun{
library( MASS )
library( VGAM, pos = 1e7, warn.conflicts = FALSE )

nObs <- 10000

# error terms (trivariate normal)
sigma <- symMatrix( c( 2, 0.7, 1.2, 1, 0.5, 1 ) )
myData <- as.data.frame( mvrnorm( nObs, c( 0, 0, 0 ), sigma ) )
names( myData ) <- c( "e0", "e1", "e2" )

# exogenous variables (indepently normal)
myData$x0 <- rnorm( nObs )
myData$x1 <- rnorm( nObs )
myData$x2 <- rnorm( nObs )

# endogenous variables
myData$y0 <-   -1.5 + 0.8 * myData$x1 + myData$e0
myData$y1 <- (  0.3 + 0.4 * myData$x1 + 0.3 * myData$x2 + myData$e1 ) > 0
myData$y2 <- ( -0.1 + 0.6 * myData$x1 + 0.7 * myData$x2 + myData$e2 ) > 0

# bivariate probit (using rhobit transformation)
bProbit <- vglm( cbind( y1, y2 ) ~ x1 + x2, family = binom2.rho,
   data = myData )
summary( bProbit )

# bivariate probit (NOT using rhobit transformation)
bProbit2 <- vglm( cbind( y1, y2 ) ~ x1 + x2, family = binom2.rho(
   lrho = "identity" ), data = myData )
summary( bProbit2 )

# inverse Mills Ratios
imr  <- invMillsRatio( bProbit )
imr2 <- invMillsRatio( bProbit2 )
all.equal( imr, imr2, tolerance = .Machine$double.eps ^ 0.25)

# tests
# E[ e0 | y1* > 0 & y2* > 0 ]
mean( myData$e0[ myData$y1 & myData$y2 ] )
mean( sigma[1,2] * imr$IMR11a + sigma[1,3] * imr$IMR11b, na.rm = TRUE )
# E[ e0 | y1* > 0 & y2* <= 0 ]
mean( myData$e0[ myData$y1 & !myData$y2 ] )
mean( sigma[1,2] * imr$IMR10a + sigma[1,3] * imr$IMR10b, na.rm = TRUE )
# E[ e0 | y1* <= 0 & y2* > 0 ]
mean( myData$e0[ !myData$y1 & myData$y2 ] )
mean( sigma[1,2] * imr$IMR01a + sigma[1,3] * imr$IMR01b, na.rm = TRUE )
# E[ e0 | y1* <= 0 & y2* <= 0 ]
mean( myData$e0[ !myData$y1 & !myData$y2 ] )
mean( sigma[1,2] * imr$IMR00a + sigma[1,3] * imr$IMR00b, na.rm = TRUE )
# E[ e0 | y1* > 0 ]
mean( myData$e0[ myData$y1 ] )
mean( sigma[1,2] * imr$IMR1X, na.rm = TRUE )
# E[ e0 | y1* <= 0 ]
mean( myData$e0[ !myData$y1 ] )
mean( sigma[1,2] * imr$IMR0X, na.rm = TRUE )
# E[ e0 | y2* > 0 ]
mean( myData$e0[ myData$y2 ] )
mean( sigma[1,3] * imr$IMRX1, na.rm = TRUE )
# E[ e0 | y2* <= 0 ]
mean( myData$e0[ !myData$y2 ] )
mean( sigma[1,3] * imr$IMRX0, na.rm = TRUE )

# estimation for y1* > 0 and y2* > 0
selection <- myData$y1 & myData$y2
# OLS estimation
ols11 <- lm( y0 ~ x1, data = myData, subset = selection )
summary( ols11 )
# heckman type estimation
heckit11 <- lm( y0 ~ x1 + IMR11a + IMR11b, data = cbind( myData, imr ),
   subset = selection )
summary( heckit11 )

# estimation for y1* > 0 and y2* <= 0
selection <- myData$y1 & !myData$y2
# OLS estimation
ols10 <- lm( y0 ~ x1, data = myData, subset = selection )
summary( ols10 )
# heckman type estimation
heckit10 <- lm( y0 ~ x1 + IMR10a + IMR10b, data = cbind( myData, imr ),
   subset = selection )
summary( heckit10 )

# estimation for y1* <= 0 and y2* > 0
selection <- !myData$y1 & myData$y2
# OLS estimation
ols01 <- lm( y0 ~ x1, data = myData, subset = selection )
summary( ols01 )
# heckman type estimation
heckit01 <- lm( y0 ~ x1 + IMR01a + IMR01b, data = cbind( myData, imr ),
   subset = selection )
summary( heckit01 )

# estimation for y1* <= 0 and y2* <= 0
selection <- !myData$y1 & !myData$y2
# OLS estimation
ols00 <- lm( y0 ~ x1, data = myData, subset = selection )
summary( ols00 )
# heckman type estimation
heckit00 <- lm( y0 ~ x1 + IMR00a + IMR00b, data = cbind( myData, imr ),
   subset = selection )
summary( heckit00 )

# estimation for y1* > 0
selection <- myData$y1
# OLS estimation
ols1X <- lm( y0 ~ x1, data = myData, subset = selection )
summary( ols1X )
# heckman type estimation
heckit1X <- lm( y0 ~ x1 + IMR1X, data = cbind( myData, imr ),
   subset = selection )
summary( heckit1X )

# estimation for y1* <= 0
selection <- !myData$y1
# OLS estimation
ols0X <- lm( y0 ~ x1, data = myData, subset = selection )
summary( ols0X )
# heckman type estimation
heckit0X <- lm( y0 ~ x1 + IMR0X, data = cbind( myData, imr ),
   subset = selection )
summary( heckit0X )

# estimation for y2* > 0
selection <- myData$y2
# OLS estimation
olsX1 <- lm( y0 ~ x1, data = myData, subset = selection )
summary( olsX1 )
# heckman type estimation
heckitX1 <- lm( y0 ~ x1 + IMRX1, data = cbind( myData, imr ),
   subset = selection )
summary( heckitX1 )

# estimation for y2* <= 0
selection <- !myData$y2
# OLS estimation
olsX0 <- lm( y0 ~ x1, data = myData, subset = selection )
summary( olsX0 )
# heckman type estimation
heckitX0 <- lm( y0 ~ x1 + IMRX0, data = cbind( myData, imr ),
   subset = selection )
summary( heckitX0 )
}
}

\keyword{models}
