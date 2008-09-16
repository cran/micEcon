summary.translogEst <- function( object, ... ) {

   object$coefTable <- coefTable( coef( object ),
         diag( vcov( object ) )^0.5, df.residual( object$lm ) )

   class( object ) <- "summary.translogEst"
   return( object )
}
