maxBHHH <- function(fn, grad=NULL, hess=NULL,
                    theta,
                    print.level=0,
                    iterlim=100,
                    ...) {
   ## HessLik   Hessian, not used, for compatibility reasons
   gradVal <- NULL
   # Save the value of gradient and use it later for hessian
   # Hessian must be called with the same parameter as gradient
   gradient <- function(theta, ...) {
      if(!is.null(grad)) {
         g <- grad(theta, ...)
      } else {
         g <- t(numericGradient(fn, theta, ...))
      }
      assign("gradVal", g, inherits=TRUE)
      return( g )
   }
   hess <- function(theta, ...) {
      g <- gradVal
      return( -t(g) %*% g )
   }
   a <- maxNR(fn, grad=gradient, hess=hess, theta=theta, iterlim=iterlim,
              print.level=print.level, ...)
   a$type = "BHHH maximisation"
   invisible(a)
}
