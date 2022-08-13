########################################################################
## This function is adapted/modified based on the plot.cv function
## from the `glmnet` package:
##
## Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
## Regularization Paths for Generalized Linear Models via Coordinate Descent.
##   Journal of Statistical Software, 33(1), 1-22.
##   URL: https://www.jstatsoft.org/v33/i01/.

##' Plot the cross-validation curve produced by cv.gcdnet
##'
##' Plots the cross-validation curve, and upper and lower standard deviation
##' curves, as a function of the \code{lambda} values used. This function is
##' modified based on the \code{plot.cv} function from the \code{glmnet}
##' package.
##'
##' A plot is produced.
##'
##' @param x fitted \code{\link{cv.gcdnet}} object
##' @param sign.lambda either plot against \code{log(lambda)} (default) or its
##'   negative if \code{sign.lambda=-1}.
##' @param \dots other graphical parameters to plot
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'
##' Maintainer: Yi Yang <yi.yang6@mcgill.ca>
##' @seealso \code{\link{cv.gcdnet}}.
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
##' # fit an elastic net penalized logistic regression with lambda2 = 1 for the
##' # L2 penalty. Use the logistic loss as the cross validation prediction loss.
##' # Use five-fold CV to choose the optimal lambda for the L1 penalty.
##' data(FHT)
##' set.seed(2011)
##' cv=cv.gcdnet(FHT$x, FHT$y, method ="logit", lambda2 = 1,
##'              pred.loss="loss", nfolds=5)
##' plot(cv)
##'
##' @importFrom graphics points abline
##' @export
plot.cv.gcdnet <- function(x, sign.lambda = 1, ...) {
  cvobj <- x
  xlab <- "log(Lambda)"
  if (sign.lambda < 0)
    xlab <- paste("-", xlab, sep = "")
  plot.args <- list(x = sign.lambda * log(cvobj$lambda), y = cvobj$cvm,
                    ylim = range(cvobj$cvupper, cvobj$cvlo), xlab = xlab,
                    ylab = cvobj$name, type = "n")
  new.args <- list(...)
  if (length(new.args))
    plot.args[names(new.args)] <- new.args
  do.call("plot", plot.args)
  error.bars(sign.lambda * log(cvobj$lambda), cvobj$cvupper,
             cvobj$cvlo, width = 0.01, col = "darkgrey")
  points(sign.lambda * log(cvobj$lambda), cvobj$cvm, pch = 20,
         col = "red")
  axis(side = 3, at = sign.lambda * log(cvobj$lambda), labels = paste(cvobj$nz),
       tick = FALSE, line = 0)
  abline(v = sign.lambda * log(cvobj$lambda.min), lty = 3)
  abline(v = sign.lambda * log(cvobj$lambda.1se), lty = 3)
  invisible()
}
