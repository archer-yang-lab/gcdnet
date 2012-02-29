cvcompute = function(mat, foldid, nlams) {
    ###Computes the weighted mean and SD within folds, and hence
    #   the se of the mean
    nfolds = max(foldid)
    outmat = matrix(NA, nfolds, ncol(mat))
    good = matrix(0, nfolds, ncol(mat))
    mat[is.infinite(mat)] = NA
    for (i in seq(nfolds)) {
        mati = mat[foldid == i, ]
        outmat[i, ] = apply(mati, 2, mean, na.rm = TRUE)
        good[i, seq(nlams[i])] = 1
    }
    N = apply(good, 2, sum)
    list(cvraw = outmat, N = N)
} 
