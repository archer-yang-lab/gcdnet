##' Print a gcdnet object
##'
##' Print a summary of the gcdnet path at each step along the path.
##'
##' The call that produced the \code{\link{gcdnet}} object is printed, followed
##' by a two-column matrix with columns \code{Df} and \code{Lambda}. The
##' \code{Df} column is the number of nonzero coefficients.
##'
##' @param x fitted \code{\link{gcdnet}} object
##' @param digits significant digits in printout
##' @param \dots additional print arguments
##' @return a two-column matrix, the first columns is the number of nonzero
##'   coefficients and the second column is \code{Lambda}.
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'
##' Maintainer: Yi Yang <yi.yang6@mcgill.ca>
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
##' m1 <- gcdnet(x = FHT$x, y = FHT$y, delta = 1, lambda2 = 0.1)
##' print(m1)
##'
##' @method print gcdnet
##' @export
print.gcdnet <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall: ", deparse(x$call), "\n\n")
  print(cbind(Df = x$df, Lambda = signif(x$lambda, digits)))
}
