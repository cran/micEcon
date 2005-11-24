linearPredictors <- function( x, ... ) {
    UseMethod("linearPredictors")
}

linearPredictors.probit <- function( x, ... ) {
   model.matrix(x) %*% x$estimate
}
