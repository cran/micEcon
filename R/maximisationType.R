maximisationType <- function(x)
    UseMethod("maximisationType")

maximisationType.default <- function(x)
    x$maximisationType

maximisationType.maximisation <- function(x)
    x$type

maximisationType.MLEstimate <- function(x)
    maximisationType(x$maxLik)
