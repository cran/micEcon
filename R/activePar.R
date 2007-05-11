## activePar: returns parameters which are free under maximisation (not fixed as constants)

activePar <- function(x, ...)
    UseMethod("activePar")

activePar.default <- function(x, ...)
    x$activePar
