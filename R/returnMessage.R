
returnMessage <- function(x, ...)
    UseMethod("returnMessage")

returnMessage.default <- function(x, ...)
    x$returnMessage

returnMessage.maximisation <- function(x, ...)
    x$message

returnMessage.maxLik <- function(x, ...)
    returnMessage(x$maxLik, ...)
