### Methods for accessing loglik value maximum likelihood estimates

loglikValue <- function(x, ...)
    UseMethod("loglikValue")

loglikValue.default <- function(x, ...)
    x$loglik

loglikValue.maxLik <- function(x, ...)
    x$maximum
