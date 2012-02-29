err = function(n, maxit, pmax) {
    if (n == 0) 
        msg = ""
    if (n > 0) {
        if (n < 7777) 
            msg = "Memory allocation error"
        if (n == 7777) 
            msg = "All used predictors have zero variance"
        if (n == 10000) 
            msg = "All penalty factors are <= 0"
        n = 1
        msg = paste("in gcdnet fortran code -", msg)
    }
    if (n < 0) {
        if (n > -10000) 
            msg = paste("Convergence for ", -n, "th lambda value not reached after maxit=", 
                maxit, " iterations; solutions for larger lambdas returned", 
                sep = "")
        if (n < -10000) 
            msg = paste("Number of nonzero coefficients along the path exceeds pmax=", 
                pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned", 
                sep = "")
        n = -1
        msg = paste("from gcdnet fortran code -", msg)
    }
    list(n = n, msg = msg)
} 
