model.frame.probit <- function (formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 
        0)]
    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        fcall$method <- "model.frame"
        fcall[[1]] <- as.name("probit")
        fcall[names(nargs)] <- nargs
        env <- environment(formula$terms)
        if (is.null(env)) 
            env <- parent.frame()
        eval(fcall, env, parent.frame())
    }
    else formula$model
}
