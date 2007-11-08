maxSANN <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=10000,
                    tol=1e-8, reltol=tol,
                    temp=10, tmax=10, parscale=rep(1, length=length(start)),
                    ...) {
   ## ... : further arguments to fn()
   ##
   ## Note: grad and hess are for compatibility only, SANN uses only fn values
   message <- function(c) {
      switch(as.character(c),
             "0" = "successful convergence",
             "10" = "degeneracy in Nelder-Mead simplex",
             "51" = "warning from the 'L-BFGS-B' method; see the corresponding component 'message' for details",
             "52" = "error from the 'L-BFGS-B' method; see the corresponding component 'message' for details"
             )
   }
   func <- function(theta, ...) {
      sum(fn(theta, ...))
   }
   gradient <- function(theta, ...) {
      ## this gradient function is used only for returning the gradient at the optimum
      if(!is.null(grad)) {
         g <- grad(theta, ...)
         if(!is.null(dim(g))) {
            if(ncol(g) > 1) {
               return(colSums(g))
            }
         } else {
            return(g)
         }
      }
      g <- numericGradient(fn, theta, ...)
      if(!is.null(dim(g))) {
         return(colSums(g))
      } else {
         return(g)
      }
   }
   hessian <- function(theta, ...) {
      ## just used for computing the final hessian, eventually using the supplied analytic information
      if(!is.null(hess)) {
         return(as.matrix(hess(theta, ...)))
      }
      return(numericHessian(fn, gradient, theta, ...))
   }
   type <- "SANN maximisation"
   parscale <- rep(parscale, length.out=length(start))
   control <- list(trace=print.level,
                    REPORT=1,
                   fnscale=-1,
                   reltol=reltol,
                    maxit=iterlim,
                   parscale=parscale,
                   temp=temp,
                   tmax=tmax)
   a <- optim(start, func, control=control, method="SANN", hessian=FALSE, ...)
   result <- list(
                   maximum=a$value,
                   estimate=a$par,
                   gradient=gradient(a$par),
                   hessian=hessian(a$par),
                   code=a$convergence,
                   message=paste(message(a$convergence), a$message),
                   last.step=NULL,
                   iterations=a$counts[1],
                   type=type)
   class(result) <- "maximisation"
   invisible(result)
}

