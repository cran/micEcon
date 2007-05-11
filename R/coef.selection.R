coef.selection <- function( object, part="full", ... ) {

   if( !( part %in% c( "full", "outcome" ) ) ) {
      stop( "argument 'part' must be either 'full' or 'outcome'" )
   }

   addToCoefNames <- function( prefix, index ) {
      if( !is.null( index ) ){
         names( coefValues )[ index ] <-
            paste( prefix, names( coefValues )[ index ], sep = "" )
      }
      return( coefValues )
   }
   if("maxLik" %in% class(object))
      coefValues <- coef.maxLik(object)
   else
       coefValues <- object$coefficients
   if( object$tobitType == 2) {
      coefValues <- addToCoefNames( "S:", object$param$index$betaS )
      coefValues <- addToCoefNames( "O:", object$param$index$betaO )
   } else if( object$tobitType == 5) {
      coefValues <- addToCoefNames( "S:",  object$param$index$betaS )
      coefValues <- addToCoefNames( "O1:", object$param$index$betaO1 )
      coefValues <- addToCoefNames( "O2:", object$param$index$betaO2 )
   }

   if( part == "outcome" ) {
      coefValues <- coefValues[ object$param$index$outcome ]
   }

   return( coefValues )
}
