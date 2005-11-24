coefTable <- function( coef, stdErr, df ) {
   result <- cbind(
      coef,
      stdErr,
      coef / stdErr,
      2 * pt( abs( coef / stdErr ), df, lower.tail = FALSE )
      )
   colnames( result ) <- c( "Estimate", "Std. Error", "t value", "Pr(>|t|)" )
   if( !is.null( names( coef ) ) ) {
      rownames( result ) <- names( coef )
   }
   return( result )
}