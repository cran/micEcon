predict.quadFuncEst <- function( object, newdata = NULL, ... ) {
   
   if( class( object$est )[1] != "lm" ) {
      warning( "predictions with panel data models may not return",
        " what you expect as it may ignore individual and/or time effects" )
   }
   
   if( is.null( newdata ) ) {
      result <- predict( object$est )
   } else {
      result <- quadFuncCalc( xNames = object$xNames, data = newdata,
         coef = coef( object ), shifterNames = object$shifterNames,
         homWeights = object$homWeights )
   }
   return( result )
}
