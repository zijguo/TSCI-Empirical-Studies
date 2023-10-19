tryCatch_E <- function(expr, ret.obj){
    # error handler
    errhandler <- function(e) {
        err <<- conditionMessage(e)
        ret.obj # Return argument ret.obj if an error occcurs.
    }
    
    # evaluate the expression
    err <- NULL
    value <- withCallingHandlers(tryCatch(expr, error = errhandler))
    
    # return a list with value, error and warning
    list(value = value, error = err)
}