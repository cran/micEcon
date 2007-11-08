
compareDerivatives <- function(f, grad, hess=NULL, t0, eps=1e-6, ...) {
### t0 - initial parameter vector
  cat("-------- compare derivatives -------- \n")
  a <- f(t0, ...)
  analytic <- grad(t0, ...)
  if(any(is.na(analytic)))
      stop("NA in analytic gradient!\n")
  if(is.null(dim(analytic))) {
    cat("Note: analytic gradient is vector.  Transforming into a matrix form\n")
    if(length(a) > 1)
        analytic <- matrix(analytic, length(analytic), 1)
                                        # Note: we assume t0 is a simple vector -> hence gradient
                                        # will be a column vector
    else
        analytic <- matrix(analytic, 1, length(analytic))
                                        # f returns a scalar -> we have row vector along t0
  }
  cat("Function value:", a, "\n")
  cat("Dim of analytic gradient:", dim(analytic), "\n")
  numeric <- numericGradient(f, t0, eps, ...)
  cat("       numeric          :", dim(numeric), "\n")
  rDiff <- (analytic - numeric)/analytic
  if(ncol(analytic) < 2) {
      a <- cbind(t0, analytic, numeric, rDiff)
      dimnames(a) <- list(param=names(t0), c("theta 0", "analytic", "numeric", "rel.diff"))
      print(a)
  }
  else {
      cat("t0\n")
      print(t0)
      cat("analytic gradient\n")
      print(analytic)
      cat("numeric gradient\n")
      print(numeric)
      cat("(anal - num)/anal\n")
      print(rDiff)
  }
  cat("Max relative difference:", max(abs(rDiff), na.rm=TRUE), "\n")
  if(!is.null(hess)) {
      cat("Comparing hessians: relative dfference\n")
    analytic <- hess(t0, ...)
    numeric <- numericGradient(grad, t0, eps, ...)
    print((analytic - numeric)/analytic)
  }
  cat("-------- END of compare derivatives -------- \n")
}
