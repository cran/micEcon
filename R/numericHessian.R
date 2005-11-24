numericHessian <- function(f, grad=NULL, t0, eps=1e-6, ...) {
   a <- f(t0, ...)
   if(is.null(grad)) {
      numericNHessian(f, t0, eps, ...)
                                        # gradient not provided -> everything numerically
   } else {
      numericGradient(grad, t0, eps, ...)
                                        # gradient is provided -> Hessian is grad grad
   }
}

numericNHessian <- function( f, t0, eps=1e-6, ...) {
   ### Numeric Hessian without gradient
   f00 <- f( t0, ...)
   eps2 <- eps*eps
   N <- length( t0)
   H <- matrix(NA, N, N)
   for( i in 1:N) {
      for( j in 1:N) {
         t01 <- t0
         t10 <- t0
         t11 <- t0
                                          # initial point
         t01[i] <- t01[i] + eps
         t10[j] <- t10[j] + eps
         t11[i] <- t11[i] + eps
         t11[j] <- t11[j] + eps
         f01 <- f( t01, ...)
         f10 <- f( t10, ...)
         f11 <- f( t11, ...)
         H[i,j] <- ( f11 - f01 - f10 + f00)/eps2
      }
   }
   return( H )
}

