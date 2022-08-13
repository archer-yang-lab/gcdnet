######################################################################
## This function is adapted/modified based on the plot function
## from the `glmnet` package:
##
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate Descent.
##   Journal of Statistical Software, 33(1), 1-22.
##   URL: https://www.jstatsoft.org/v33/i01/.

##' Plot coefficients from a "gcdnet" object
##'
##' Produces a coefficient profile plot of the coefficient paths for a fitted
##' \code{\link{gcdnet}} object. This function is modified based on the
##' \code{plot} function from the \code{glmnet} package.
##'
##' A coefficient profile plot is produced.
##'
##' @param x fitted \code{\link{gcdnet}} model
##' @param xvar what is on the X-axis. \code{"norm"} plots against the L1-norm
##'   of the coefficients, \code{"lambda"} against the log-lambda sequence.
##' @param color if \code{TRUE}, plot the curves with rainbow colors.
##'   \code{FALSE} is gray colors. Default is \code{FALSE}
##' @param label if \code{TRUE}, label the curves with variable sequence
##'   numbers. Default is \code{FALSE}
##' @param \dots other graphical parameters to plot
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr Maintainer: Yi Yang
##'   <yi.yang6@mcgill.ca>
##' @references Yang, Y. and Zou, H. (2012).
##'   "An Efficient Algorithm for Computing The HHSVM and Its Generalizations."
##'   \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
##'   BugReport: \url{https://github.com/emeryyi/gcdnet}\cr
##'
##'   Gu, Y., and Zou, H. (2016).
##'   "High-dimensional generalizations of asymmetric least squares regression and their applications."
##'   \emph{The Annals of Statistics}, 44(6), 2661–2694.\cr
##'
##'   Friedman, J., Hastie, T., and Tibshirani, R. (2010).
##'   "Regularization paths for generalized linear models via coordinate descent."
##'   \emph{Journal of Statistical Software, 33, 1.}\cr
##'   \url{https://www.jstatsoft.org/v33/i01/}
##' @keywords models regression
##' @examples
##'
##' data(FHT)
##' m1 <- gcdnet(x = FHT$x,y = FHT$y)
##' par(mfrow = c(1,3))
##' plot(m1) # plots against the L1-norm of the coefficients
##' plot(m1,xvar = "lambda",label = TRUE) # plots against the log-lambda sequence
##' plot(m1,color = TRUE)
##'
##' @importFrom graphics matplot axis text
##' @importFrom grDevices gray.colors rainbow
##' @importFrom stats approx
##' @export
plot.gcdnet <- function(x, xvar = c("norm", "lambda"), color = FALSE,
                        label = FALSE, ...) {
  beta <- x$beta
  lambda <- x$lambda
  df <- x$df
  xvar <- match.arg(xvar)
  ## beta should be in 'dgCMatrix' format
  which <- nonzero(beta)
  beta <- as.matrix(beta[which, ])
  xvar <- match.arg(xvar)
  switch(xvar, norm = {
    index <- apply(abs(beta), 2, sum)
    iname <- "L1 Norm"
  }, lambda = {
    index <- log(lambda)
    iname <- "Log Lambda"
  })
  xlab <- iname
  ylab <- "Coefficients"
  dotlist <- list(...)
  type <- dotlist$type
  if (is.null(type)) {
    if (color == FALSE)
      matplot(index, t(beta), lty = 1, xlab = xlab,
              ylab = ylab, type = "l", pch = 500,
              col = rainbow(12, start = 0.7, end = 0.95),
              ...) else matplot(index, t(beta), lty = 1,
                                xlab = xlab, ylab = ylab,
                                type = "l", pch = 500, ...)
  } else matplot(index, t(beta), lty = 1, xlab = xlab, ylab = ylab, ...)
  atdf <- pretty(index)
  prettydf <- trunc(approx(x = index, y = df, xout = atdf,
                           rule = 2)$y)
  axis(3, at = atdf, labels = prettydf, cex.axis = 1, tcl = NA)
  if (label) {
    nnz <- length(which)
    xpos <- max(index)
    pos <- 4
    if (xvar == "lambda") {
      xpos <- min(index)
      pos <- 2
    }
    xpos <- rep(xpos, nnz)
    ypos <- beta[, ncol(beta)]
    text(xpos, ypos, paste(which), cex = 0.5, pos = pos)
  }
}
