elas.aidsEst <- function( object, method = NULL, ... ) {

   # specify default value for argument method
   if( is.null ( method ) ) {
      if( substr( object$method, 1, 2 ) == "LA" ) {
         method <- "Ch"
      } else if( substr( object$method, 1, 2 ) %in% c( "MK", "IL" ) ) {
         method <- "AIDS"
      }
   }

   # test reasonability of argument method
   if( substr( object$method, 1, 2 ) %in% c( "MK", "IL" ) &&
         method != "AIDS" ) {
      warning( paste( "It does not make sense to calculate the elasticities",
         " of a (non-linear) AIDS model with method '", method, "'",
         sep = "" ) )
   }

   # to avoid warning message in aidsElas
   if( method %in% c( "Go", "Ch", "EU" ) ) {
      object$pMeans <- NULL
   }

   # calculate demand elasticities
   result  <- aidsElas( coef = object$coef,
      shares = object$wMeans, prices = object$pMeans,
      method = method,
      priceNames = object$priceNames,
      coefVcov = object$coef$allcov, df = object$est$df, ... )

   return( result )
}

