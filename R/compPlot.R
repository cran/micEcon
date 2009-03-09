compPlot <- function( x, y, lim = NULL, ... ) {

   xyRange <- range( x, y, na.rm = TRUE, finite = TRUE )

   if( is.null( lim ) ) {
      lim <- xyRange
   } else {
      if( length( lim ) != 2 ) {
         stop( "argument 'lim' must be a vector of two elements" )
      }
      if( is.na( lim[1] ) ) {
         lim[1] <- xyRange[1]
      }
      if( is.na( lim[2] ) ) {
         lim[2] <- xyRange[2]
      }
      if( lim[1] >= lim[2] ) {
         stop( "the first element of argument 'lim' must be smaller",
            " than the second element" )
      }
      if( lim[1] > xyRange[1] |  lim[2] < xyRange[2] ) {
         warning( "some data points are outside the print area" )
      }
   }

   plot( x, y, xlim = lim, ylim = lim, ... )

   lines( 1.5 * lim - 0.5 * lim[c(2,1)],
      1.5 * lim - 0.5 * lim[c(2,1)] )

   invisible( xyRange )
}
