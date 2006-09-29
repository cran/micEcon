vcov.maxLik <- function( object, lambdatol=1e-6, ... ) {
   ## object      object of class "maxLik"
   ## lambdatol   minimum eigenvalue where the hessian is still inverted
   ## 
   ## RESULTS:
   ## variance-covariance matrix
   ## 
   NParam <- length(object$estimate)
   if(!is.null(object$acivePar)) {
      activePar <- object$activeParp
   } else {
      activePar <- rep(TRUE, NParam)
   }
   if(object$code < 100) {
      if(min(abs(eigen(object$Hessian[activePar,activePar],
            symmetric=TRUE, only.values=TRUE)$values)) > lambdatol) {
         varcovar <- matrix(0, NParam, NParam)
         varcovar[activePar,activePar] <-
             solve(-object$Hessian[activePar,activePar])
         return(varcovar)
      }
   }
   NULL
}
