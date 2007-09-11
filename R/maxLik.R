maxLik <- function(logLik, grad=NULL, hess=NULL, start,
                   method="Newton-Raphson",
                   ...) {
   ## Maximum Likelihood estimation.
   ##
   ## Newton-Raphson maximisation
   ## Parameters:
   ## logLik     log-likelihood function.  First argument must be the vector of parameters.
   ## grad       gradient of log-likelihood.  If NULL, numeric gradient is used.  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(nObs, nParam).  In this case the rows are simply
   ##                 summed (useful for maxBHHH).
   ## hess       Hessian function (numeric used if NULL)
   ## start      initial vector of parameters (eventually w/names)
   ## method     maximisation method (Newton-Raphson)
   ## ...        additional arguments for the maximisation routine
   ##
   ## RESULTS:
   ## list of class c("maxLik", "maximisation").  This is in fact equal to class "maximisation", just the
   ## methods are different.
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
   maxRoutine <- switch(method,
                        "Newton-Raphson" =,
                        "newton-raphson" =,
                        "NR" =,
                        "nr" = maxNR,
                        "BFGS" =,
                        "bfgs" = maxBFGS,
                        "BHHH" =,
                        "bhhh" = maxBHHH,
                        stop( "Maxlik: unknown maximisation method ", method )
                        )
   result <- maxRoutine(fn=logLik, grad=grad, hess=hess, start=start, ...)
   class(result) <- c("maxLik", class(result))
   result
}
