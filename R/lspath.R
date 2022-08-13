lspath <- function(x, y, nlam, flmin, ulam, isd, intr, eps, dfmax, pmax, jd, pf,
                   pf2, maxit, lam2, nobs, nvars, vnames) {
################################################################################
  ## data setup
  y <- as.double(y)
  storage.mode(x) <- "double"
################################################################################
  ## call Fortran core
  fit <- .Fortran("lslassoNET", lam2, nobs, nvars, x, as.double(y), jd, pf, pf2,
                  dfmax, pmax, nlam, flmin, ulam, eps, isd, intr, maxit,
                  nalam = integer(1), b0 = double(nlam),
                  beta = double(pmax * nlam), ibeta = integer(pmax),
                  nbeta = integer(nlam), alam = double(nlam),
                  npass = integer(1), jerr = integer(1), PACKAGE = "gcdnet")
################################################################################
  ## output
  outlist <- getoutput(fit, maxit, pmax, nvars, vnames)
  outlist <- c(outlist, list(npasses = fit$npass, jerr = fit$jerr))
  class(outlist) <- c("lspath")
  outlist
}
