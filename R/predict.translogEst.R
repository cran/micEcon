predict.translogEst <- function( object, newdata = NULL,
   dataLogged = object$dataLogged, ... ) {

   if( !is.null( newdata ) ) {
      if( dataLogged ) {
         logData <- newdata
      } else {
         logData <- logDataSet( data = newdata,
            varNames = c( object$xNames ),
            varNamesNum = object$shifterNames )
      }
   } else {
      logData <- NULL
   }
   
   result <- predict.quadFuncEst( object, newdata = logData, ... )
   
   if( !dataLogged ) {
      result <- exp( result )
   }

   return( result )
}
