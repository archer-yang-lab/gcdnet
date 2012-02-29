getoutput = function(fit, maxit, pmax, nvars, vnames) {
    nalam = fit$nalam
    nbeta = fit$nbeta[seq(nalam)]
    nbetamax = max(nbeta)
    lam = fit$alam[seq(nalam)]
    stepnames = paste("s", seq(nalam) - 1, sep = "")
    errmsg = err(fit$jerr, maxit, pmax)
    switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), 
        `-1` = print(errmsg$msg, call. = FALSE))
    dd = c(nvars, nalam)
    if (nbetamax > 0) {
        beta = matrix(fit$beta[seq(pmax * nalam)], pmax, nalam)[seq(nbetamax), 
            , drop = FALSE]
        df = apply(abs(beta) > 0, 2, sum)
        ja = fit$ibeta[seq(nbetamax)]
        oja = order(ja)
        ja = rep(ja[oja], nalam)
        ibeta = cumsum(c(1, rep(nbetamax, nalam)))
        beta = new("dgCMatrix", Dim = dd, Dimnames = list(vnames, 
            stepnames), x = as.vector(beta[oja, ]), p = as.integer(ibeta - 
            1), i = as.integer(ja - 1))
    }
    else {
        beta = zeromat(nvars, nalam, vnames, stepnames)
        df = rep(0, nalam)
    }
    b0 = fit$b0
    if (!is.null(b0)) {
        b0 = b0[seq(nalam)]
        names(b0) = stepnames
    }
    list(b0 = b0, beta = beta, df = df, dim = dd, lambda = lam)
} 
