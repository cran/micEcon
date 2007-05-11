print.summary.maximisation <- function( x, ... ) {
   summary <- x
   cat("--------------------------------------------\n")
   cat(summary$type, "\n")
   cat("Number of iterations:", summary$iterations, "\n")
   cat("Return code:", summary$code, "\n")
   cat(summary$message, "\n")
   if(!is.null(summary$unsucc.step)) {
      cat("Last (unsuccessful) step: function value", summary$unsucc.step$value,
         "\n")
      print(summary$unsucc.step$parameters)
   }
   if(!is.null(summary$estimate)) {
      cat("Function value:", summary$maximum, "\n")
      cat("Estimates:\n")
      print(summary$estimate)
      if(!is.null(summary$hessian)) {
         cat("Hessian:\n")
         print(summary$hessian)
      }
   }
   cat("--------------------------------------------\n")
}

summary.maximisation <- function(object, hessian=FALSE, unsucc.step=FALSE,
   ... ) {
   ## The object of class "maximisation" should include following components:
   ## maximum    : function value at optimum
   ## estimate   : matrix, estimated parameter values and gradient at optimum
   ## hessian    :           hessian
   ## code       : code of convergence
   ## message    : message, description of the code
   ## last.step  : information about last step, if unsuccessful
   ## iterations : number of iterations
   ## type       : type of optimisation
   ##
   nParam <- length(object$estimate)
   if(!is.null(object$acivePar)) {
      activePar <- object$activePar
   } else {
      activePar <- rep(TRUE, nParam)
   }
   if(object$code == 3 & unsucc.step) {
      a <- cbind(object$last.step$theta0, object$last.step$theta1)
      dimnames(a) <- list(parameter=object$names,
                        c("current par", "new par"))
      unsucc.step <- list(value=object$last.step$f0,
                        parameters=a)
   } else {
      unsucc.step <- NULL
   }
   estimate <- cbind("estimate"=object$estimate, "gradient"=object$gradient)
   if(hessian) {
      H <- object$hessian
   }
   else {
      H <- NULL
   }
   summary <- list(type=object$type,
                  iterations=object$iterations,
                  code=object$code,
                  message=object$message,
                   unsucc.step=unsucc.step,
                   maximum=object$value,
                  estimate=estimate,
                   hessian=H)
   class(summary) <- "summary.maximisation"
   summary
}
