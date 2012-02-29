nonzero = function(beta, bystep = FALSE) {
	ns = ncol(beta)
    ##beta should be in 'dgCMatrix' format
    if (nrow(beta) == 1) {
        if (bystep) 
            apply(beta, 2, function(x) if (abs(x) > 0) 
                1
            else NULL)
        else {
            if (any(abs(beta) > 0)) 
                1
            else NULL
        }
    }
    else {
        beta = t(beta)
        which = diff(beta@p)
        which = seq(which)[which > 0]
        if (bystep) {
            nzel = function(x, which) if (any(x)) 
                which[x]
            else NULL
			beta = abs(as.matrix(beta[, which])) > 0
			if(ns==1) apply(beta, 2, nzel, which)
            else apply(beta, 1, nzel, which)
        }
        else which
    }
} 
