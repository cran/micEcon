print.summary.probit <- function( x, ... ) {
   cat("--------------------------------------------\n")
   cat("Probit binary choice model/Maximum Likelihood estimation\n")
   cat(x$type, ", ", x$iterations, " iterations\n", sep="")
   cat("Return code ", x$code, ": ", x$message, "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$loglik, "\n")
      cat(x$nObs, " observations (", x$N0, " zeros and ", x$N1, " ones) and ",
          x$NActivePar, " free parameters (df = ",
          x$nObs - x$NActivePar, ")\n", sep="")
      cat("Estimates:\n")
      print(x$estimate)
   }
   cat("Significance test:\n")
   cat("chi2(", x$LRT$df, ") = ", x$LRT$LRT, " (p=", x$LRT$pchi2, ")\n", sep="")
   cat("--------------------------------------------\n")
}
