aidsSystem <- function( nGoods, LA = TRUE ) {
  if( LA ) {
    system <- list()
    for(i in 1:( nGoods - 1 ) ) {
      system[[ i ]] <- paste( "w", as.character( i ), " ~ lxtr", sep = "" )
      for( j in 1:nGoods ) {
        system[[ i ]] <- paste( system[[ i ]], " + lp",
           as.character( j ), sep = "" )
      }
      system[[ i ]] <- as.formula( system[[ i ]] )
    }
  }
  return( system )
}
