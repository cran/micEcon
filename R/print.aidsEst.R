print.aidsEst <- function( x, ... ) {
   cat( "\nDemand analysis with the Almost Ideal " )
   cat( "Demand System (AIDS)\n" )
   cat( "Estimation Method: " )
   cat( .aidsEstMethod( x$method, x$px ) )
   cat( "Estimated Coefficients:\n" )
   cat( "alpha\n" )
   print( x$coef$alpha )
   cat( "beta\n" )
   print( x$coef$beta )
   cat( "gamma\n" )
   print( x$coef$gamma )
   if( !is.null( x$coef$delta ) ){
      cat( "delta\n" )
      print( x$coef$delta )
   }
   invisible( x )
}
