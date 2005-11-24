summary.tobit2 <- function(object, ...) {
   ## object      object of class "tobit2"
   ## ...         additional arguments for "summary.maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.tobit2" with following components:
   ## maximum    : function value at optimum
   ## estimateProbit : probit part of the estimate
   ## estimateEquation : equation -"-
   ## estimateEps      : sigma & rho part
   ## + additional components of "summary.maxLik"
   ## 
   sl <- summary.maxLik(object, ...)
   iZ <- seq(length=object$NZ)
   iX <- tail(iZ, n=1) + seq(length=object$NX)
   iSigma <- tail(iX, n=1) + 1
   iRho <- iSigma + 1
   s <- c(sl,
          estimateProbit=list(sl$estimate[iZ,]),
          estimateEquation=list(sl$estimate[iX,]),
          estimateEps=list(sl$estimate[c(iSigma, iRho),]),
          NObs=object$NObs, N1=object$N1, N2=object$N2, NZ=object$NZ, NX=object$NX, df=object$df
          )
   class(s) <- c("summary.tobit2", class(sl))
   s
}

print.summary.tobit2 <- function(x, ...) {
   cat("--------------------------------------------\n")
   cat("Tobit 2 selection model/Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$estimate$value, "\n")
      cat(x$NObs, " observations (", x$N1, " censored and ", x$N2, " observed) and ",
          x$NActivePar, " free parameters (df =",
          x$NObs - x$NActivePar, ")\n", sep="")
      cat("\nProbit selection equation:\n")
      print(x$estimateProbit)
      cat("\nOLS equation:\n")
      print(x$estimateEquation)
      cat("\nError terms data:\n")
      print(x$estimateEps)
      if(!is.null(x$Hessian)) {
         cat("Hessian:\n")
         print(x$Hessian)
      }
   }
   cat("--------------------------------------------\n")
}
