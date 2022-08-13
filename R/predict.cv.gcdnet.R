##' Make predictions from a "cv.gcdnet" object.
##'
##' This function makes predictions from a cross-validated gcdnet model, using
##' the stored \code{"gcdnet.fit"} object, and the optimal value chosen for
##' \code{lambda}.
##'
##' This function makes it easier to use the results of cross-validation to
##' make a prediction.
##'
##' @param object fitted \code{\link{cv.gcdnet}} object.
##' @param newx matrix of new values for \code{x} at which predictions are to be
##'   made. Must be a matrix. See documentation for \code{predict.gcdnet}.
##' @param s value(s) of the penalty parameter \code{lambda} at which
##'   predictions are required. Default is the value \code{s="lambda.1se"}
##'   stored on the CV object. Alternatively \code{s="lambda.min"} can be used.
##'   If \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
##'   used.
##' @param \dots not used. Other arguments to predict.
##' @return The object returned depends the \dots{} argument which is passed on
##'   to the \code{\link{predict}} method for \code{\link{gcdnet}} objects.
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'
##' Maintainer: Yi Yang <yi.yang6@mcgill.ca>
##' @seealso \code{\link{cv.gcdnet}}, and \code{\link{coef.cv.gcdnet}} methods.
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
##' set.seed(2011)
##' cv=cv.gcdnet(FHT$x, FHT$y, lambda2 = 1, pred.loss="misclass",
##'              lambda.factor=0.05, nfolds=5)
##' pre = predict(cv$gcdnet.fit, newx = FHT$x, s = cv$lambda.1se,
##'               type = "class")
##'
##' @export
predict.cv.gcdnet <- function(object, newx, s = c("lambda.1se", "lambda.min"),
                              ...) {
  if (is.numeric(s))
    lambda <- s else if (is.character(s)) {
                  s <- match.arg(s)
                  lambda <- object[[s]]
                } else stop("Invalid form for s")
  predict(object$gcdnet.fit, newx, s = lambda, ...)
}
