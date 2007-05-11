## probit
vcov.probit <- function(object, ...) {
  result <- vcov.maxLik( object )
  if(!is.null(result))
      rownames( result ) <- colnames( result ) <- names( object$estimate )
  return( result )
}


## maxLik
vcov.maxLik <- function(object, ...) {
   ## if exists $varcovar, take it
   if(!is.null(object$varcovar))
       return(object$varcovar)
   ## otherwise invert hessian
   activePar <- activePar(object)
   if(min(abs(eigen(hessian(object)[activePar,activePar],
                    symmetric=TRUE, only.values=TRUE)$values)) > 1e-6) {
      varcovar <- matrix(0, nParam(object), nParam(object))
      varcovar[activePar,activePar] <- solve(-hessian(object)[activePar,activePar])
   }
   else
       varcovar <- NULL
   varcovar
}
