maxBFGS <- function(fn, grad=NULL, hess=NULL,
                    start,
                    print.level=0,
                    iterlim=200,
                    tol=1e-8, reltol=tol,
                    ...) {
   ## ... : further arguments to fn() and grad()
   message <- function(c) {
      switch(as.character(c),
               "0" = "successful convergence",
               "10" = "degeneracy in Nelder-Mead simplex"
               )
   }
   func <- function(theta, ...) {
      sum(fn(theta, ...))
   }
   gradient <- function(theta, ...) {
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
   type <- "BFGS maximisation"
   control <- list(trace=print.level,
                    REPORT=1,
                    fnscale=-1,
                   reltol=reltol,
                    maxit=iterlim)
   a <- optim(start, func, gr=gradient, control=control, method="BFGS",
      hessian=TRUE, ...)
   result <- list(
                   maximum=a$value,
                   estimate=a$par,
                   gradient=gradient(a$par),
                   hessian=a$hessian,
                   code=a$convergence,
                   message=paste(message(a$convergence), a$message),
                   last.step=NULL,
                   iterations=a$counts,
                   type=type)
   class(result) <- "maximisation"
   invisible(result)
}

