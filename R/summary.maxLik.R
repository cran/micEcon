
coef.maxLik <- function(object, ...)
    object$estimate

print.summary.maxLik <- function( x, ... ) {
   cat("--------------------------------------------\n")
   cat("Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$loglik, "\n")
      cat(x$NActivePar, " free parameters\n")
      cat("Estimates:\n")
      print(x$estimate)
      if(!is.null(x$Hessian)) {
         cat("Hessian:\n")
         print(x$Hessian)
      }
   }
   cat("--------------------------------------------\n")
}

summary.maxLik <- function( object, hessian=FALSE, ... ) {
   ## object      object of class "maxLik"
   ## hessian     logical, whether to include hessian in summary
   ## 
   ## RESULTS:
   ## list of class "summary.maxLik" with following components:
   ## maximum    : function value at optimum
   ## estimate   : estimated parameter values at optimum
   ## gradient   :           gradient at optimum
   ## hessian    :           hessian
   ## code       : code of convergence
   ## message    : message, description of the code
   ## iterations : number of iterations
   ## type       : type of optimisation
   ##
   result <- object$maximisation
   NParam <- length(coef <- coefficients(object))
   if(!is.null(object$activePar)) {
      activePar <- object$activePar
   } else {
      activePar <- rep(TRUE, NParam)
   }
   if(object$code < 100) {
      if(min(abs(eigen(object$Hessian[activePar,activePar],
                       symmetric=TRUE, only.values=TRUE)$values)) > 1e-6) {
         varcovar <- matrix(0, NParam, NParam)
         varcovar[activePar,activePar] <-
             solve(-object$Hessian[activePar,activePar])
         hdiag <- diag(varcovar)
         if(any(hdiag < 0)) {
            warning("Diagonal of variance-covariance matrix not positive!\n")
         }
         stdd <- sqrt(hdiag)
         t <- coef/stdd
         p <- 2*pnorm( -abs( t))
      } else {
         stdd <- 0
         t <- 0
         p <- 0
      }
      results <- cbind(coef, stdd, t, "P(|b| > t)"=p)
      Hess <- NULL
      if(hessian) {
         Hess <- object$Hessian
      }
   } else {
      results <- NULL
      Hess <- NULL
   }
   summary <- list(type=object$type,
                   iterations=object$iterations,
                   code=object$code,
                   message=object$message,
                   loglik=object$maximum,
                   estimate=results,
                   Hessian=Hess,
                   activePar=object$activePar,
                   NActivePar=sum(object$activePar))
   class(summary) <- "summary.maxLik"
   summary
}
