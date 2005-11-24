numericGradient <- function(f, t0, eps=1e-6, ...) {
   ### numeric gradient of a vector-valued function
   ### f    : function, return Nval x 1 vector of values
   ### t0   : Npar x 1 vector of parameters
   ### return:
   ### NparxNval matrix, gradient
   Npar <- length(t0)
   Nf <- length(f0 <- f(t0, ...))
   grad <- matrix(NA, Npar, Nf)
   colnames(grad) <- names(f0)
   row.names(grad) <- names(t0)
   for(i in 1:Npar) {
      t2 <- t1 <- t0
      t1[i] <- t0[i] - eps/2
      t2[i] <- t0[i] + eps/2
      grad[i,] <- (f(t2, ...) - f(t1, ...))/eps
   }
   return( grad )
}

