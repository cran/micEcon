maxBHHH <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=100,
                    ...) {
   ## hess:   Hessian, not used, for compatibility with the other methods
   gradVal <- NULL
   # Save the value of gradient and use it later for hessian
   # Hessian must be called with the same parameter as gradient
   gradient <- function(theta, ...) {
      if(!is.null(grad)) {
         g <- grad(theta, ...)
      } else {
         g <- numericGradient(fn, theta, ...)
      }
      ## Ensure g is suitable for information equality approximation
      if(is.null(dim(g)))
          g <- matrix(g)
      if(!(dim(g)[1] > length(theta) & dim(g)[2] == length(theta))) {
         stop(paste("Gradient matrix must have at least as many rows and exactly as many columns as the number of parameters.\n",
                    "Currently", length(theta), "parameters but the gradient is", dim(g)[1], "x", dim(g)[2]))
      }
      ##
      assign("gradVal", g, inherits=TRUE)
      return( g )
   }
   hess <- function(theta, ...) {
      g <- gradVal
      return( -t(g) %*% g )
   }
   a <- maxNR(fn, grad=gradient, hess=hess, start=start, iterlim=iterlim,
              print.level=print.level, ...)
   a$type = "BHHH maximisation"
   invisible(a)
}
