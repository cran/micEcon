
maxNR <- function(fn, grad=NULL, hess=NULL, start,
                  print.level=0,
                  tol=1e-6, gradtol=1e-6, steptol=1e-10,
                  lambdatol=1e-6,
                  qrtol=1e-10,
                  iterlim=15,
                  constPar=NULL,
                  activePar=rep(TRUE, nParam),
                  ...) {
   ## Newton-Raphson maximisation
   ## Parameters:
   ## fn          - the function to be minimized.  Returns either scalar or
   ##               vector value with possible attributes constPar and
   ##               constVal
   ## grad        - gradient function (numeric used if missing).  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(M, nParam), where M is arbitrary.  In this case the
   ##                 rows are simply summed (useful for maxBHHH).
   ## hess        - hessian function (numeric used if missing)
   ## start       - initial parameter vector (eventually w/names)
   ## steptol     - minimum step size
   ## lambdatol   - max lowest eigenvalue when forcing pos. definite H
   ## qrtol       - tolerance for qr decomposition
   ## ...         - extra arguments for fn()
   ## The stopping criteria
   ## tol         - maximum allowed difference between sequential values
   ## gradtol     - maximum allowed norm of gradient vector
   ## iterlim     - maximum # of iterations
   ## constPar    - NULL or an index vector -- which parameters are taken as
   ##               constants
   ##
   ## RESULTS:
   ## a list of class "maximisation":
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success:
   ##             1 - gradient close to zero
   ##             2 - successive values within tolerance limit
   ##             3 - could not find a higher point (step error)
   ##             4 - iteration limit exceeded
   ##             100 - initial value out of range
   ## message     character message describing the code
   ## last.step   only present if code == 3 (step error).  A list with following components:
   ##             teeta0 - parameetrid, millel viga tuli
   ##             f0 - funktsiooni väärtus nende parameetritega (koos
   ##                  gradiendi ja hessi maatriksiga)
   ##             teeta1 - ruutpolünoomi järgi õige uus parameetri väärtus
   ##             activePar - logical vector, which parameters are active (not constant)
   ## activePar   logical vector, which parameters were treated as free (resp fixed)
   ## iterations  number of iterations
   ## type        "Newton-Raphson maximisation"
   maximisation.message <- function(code) {
      message <- switch(code,
         "1" = "gradient close to zero. May be a solution",
         "2" = paste("successive function values within tolerance",
                     "limit.\n May be a solution"),
         "3" = paste("Last step could not find a value above the",
                     "current.\nMay be near a solution"),
         "4" = "Iteration limit exceeded.",
         "100" = "Initial value out of range.",
         paste("Code", code))
      return(message)
   }
   max.eigen <- function( M) {
      ## return maximal eigenvalue of (symmetric) matrix
      val <- eigen(M, symmetric=TRUE, only.values=TRUE)$values
      val[1]
      ## L - eigenvalues in decreasing order, [1] - biggest in abs value
   }
   func <- function(theta, ...) {
      f <- fn(theta, ...)
      sf <- sum(f)
      mostattributes(sf) <- attributes(f)
      sf
   }
   gradient <- function(theta, ...) {
      if(!is.null(grad)) {  # use user-supplied if present
         gr <- grad(theta, ...)
      } else {
         gr <- numericGradient(fn, theta, ...)
                                        # Note we need nObs rows x nParam cols
      }
      ## Now check if the gradient is vector or matrix...
      if(!is.null(dim(gr))) {
         return(colSums(gr))
      } else {
         ## ... or vector if only one parameter
         if(length(gr) > nParam) {
            return(sum(gr))
         }
      }
      return(gr)
   }
   hessian <- function(theta, ...) {
      if(!is.null(hess)) {
         return(as.matrix(hess(theta, ...)))
      }
      return(numericHessian(fn, gradient, theta, ...))
   }
   ## -------------------------------------------------
   maximisation.type <- "Newton-Raphson maximisation"
   nimed <- names(start)
   nParam <- length(start)
   I <- diag(rep(1, nParam))     # I is unit matrix
   activePar[constPar] <- FALSE
   start1 <- start
   iter <- 0
   f1 <- func(start1, ...)
   if(print.level > 2) {
      cat("Initial function value:", f1, "\n")
   }
   if(is.na( f1) | is.infinite(f1)) {
      result <- list(code=100, message=maximisation.message("100"),
                     iterations=0,
                     type=maximisation.type)
      class(result) <- "maximisation"
      return(result)
   }
   G1 <- gradient(start, ...)
   if(print.level > 2) {
      cat("Initial gradient value:\n")
      print(G1)
   }
   if(any(is.na(G1))) {
      stop("Na in the initial gradient")
   }
   if(length(G1) != nParam) {
      stop( "length of gradient (", length(G1),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   H1 <- hessian(start)
   if(any(is.na(H1))) {
      stop("NA in the initial Hessian")
   }
   if( print.level > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(f1), "\n")
      a <- cbind(start, G1, as.integer(activePar))
      dimnames(a) <- list(nimed, c("parameter", "initial gradient",
                                          "free"))
      print(a)
      cat( "Condition number of the hessian:", kappa( H1), "\n")
      if( print.level > 3) {
         print( H1)
      }
   }
   repeat {
      iter <- iter + 1
      lambda <- 0
      start0 <- start1
      f0 <- f1
      G0 <- G1
      if(any(is.na(G0))) {
         stop("NA in gradient (at the iteration start)")
      }
      H0 <- H1
      if(any(is.na(H0))) {
         stop("NA in Hessian (at the iteration start)")
      }
      step <- 1
      H <- H0
      ## check whether hessian is positive definite
      while((me <- max.eigen( H[activePar,activePar,drop=FALSE])) >= -lambdatol |
         (qRank <- qr(H[activePar,activePar], tol=qrtol)$rank) < sum(activePar)) {
                                        # maximum eigenvalue -> negative definite
                                        # qr()$rank -> singularity
         H <- H - (abs(me) + lambdatol)*I
                                        # how to make it better?
      }
      amount <- vector("numeric", nParam)
      amount[activePar] <- qr.solve(H[activePar,activePar,drop=FALSE],
                                    G0[activePar], tol=qrtol)
      start1 <- start0 - step*amount
      f1 <- func(start1, ...)
      ## Find out the constant parameters
      constPar <- attr(f1, "constPar")
      if(!is.null(constPar)) {
         if(any(is.na(constPar))) {
            stop("NA in the list of constants")
         }
         activePar <- rep(TRUE, nParam)
         activePar[constPar] <- FALSE
         if(!is.null(attr(f1, "constVal"))) {
            start1[constPar] <- attr(f1, "constVal")
               # put constants into start.  func() should
               # already use them
         }
      }
      if(is.null(newVal <- attr(f1, "newVal"))) {
         while( is.na( f1) || ( ( f1 < f0) && ( step > steptol))) {
                                        # We end up in a NA or a higher value.
                                        # try smaller step
            step <- step/2
            start1 <- start0 - step*amount
            f1 <- func(start1, ...)
            ## Find out the constant parameters -- these may be other than
            ## with full step
            constPar <- attr(f1, "constPar")
            if(!is.null(constPar)) {
               if(any(is.na(constPar))) {
                  stop("NA in the list of constants")
               }
               activePar[constPar] <- FALSE
               if(!is.null(attr(f1, "constVal"))) {
                  start1[constPar] <- attr(f1, "constVal")
               }
            }
         }
      } else {
         start1[newVal$index] <- newVal$val
      }
      G1 <- gradient(start1, ...)
      if(any(is.na(G1))) {
         cat("Iteration", iter, "\n")
         cat("Parameter:\n")
         print(start1)
         stop("NA in gradient")
      }
      H1 <- hessian(start1, ...)
      if( print.level > 1) {
        cat( "-----Iteration", iter, "-----\n")
      }
      if(print.level > 2) {
         cat( "lambda ", lambda, " step", step, " fcn value:",
            formatC(as.vector(f1), digits=8, format="f"),  "\n")
         a <- cbind(amount, start1, G1, as.integer(activePar))
         dimnames(a) <- list(names(start0), c("amount", "new param",
                                             "new gradient", "free"))
         print(a)
         cat( "Condition number of the hessian:",
            kappa(H1[activePar,activePar,drop=FALSE]), "\n")
         if( print.level > 3) {
            print( H1)
         }
      }
      if( step < steptol) {
         code <- 3; break
      }
      if( iter > iterlim) {
         code <- 4; break
      }
      if( sqrt( t(G1[activePar])%*%G1[activePar]) < gradtol) {
         code <-1; break
      }
      if(is.null(newVal) & f1 - f0 < tol) {
         code <- 2; break
      }
   }
   if( print.level > 0) {
      cat( "--------------\n")
      cat( maximisation.message( code), "\n")
      cat( iter, " iterations\n")
      cat( "estimate:", start1, "\n")
      cat( "Function value:", f1, "\n")
   }
   if( code == 3) {
      samm <- list(theta0=start0, f0=f0, start1=start1)
   } else {
      samm <- NULL
   }
   names(start1) <- nimed
   result <-list(
                  maximum=as.vector( f1),
                  estimate=start1,
                  gradient=G1,
                 hessian=H1,
                  code=code,
                  message=maximisation.message( code),
                  last.step=samm,
                                        # only when could not find a
                                        # lower point
                  activePar=activePar,
                  iterations=iter,
                  type=maximisation.type)
   class(result) <- c("maximisation", class(result))
   invisible(result)
}

returnCode.maximisation <- function(x, ...)
    x$code
