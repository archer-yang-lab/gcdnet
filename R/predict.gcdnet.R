##' Make predictions from a "gcdnet" object
##'
##' Similar to other predict methods, this functions predicts fitted values and
##' class labels from a fitted \code{\link{gcdnet}} object.
##'
##' \code{s} is the new vector at which predictions are requested. If \code{s}
##' is not in the lambda sequence used for fitting the model, the
##' \code{predict} function will use linear interpolation to make predictions.
##' The new values are interpolated using a fraction of predicted values from
##' both left and right \code{lambda} indices.
##'
##' @aliases predict.gcdnet predict.hsvmpath predict.sqsvmpath predict.logitpath
##'   predict.lspath predict.erpath
##'
##' @param object fitted \code{\link{gcdnet}} model object.
##' @param newx matrix of new values for \code{x} at which predictions are to be
##'   made. NOTE: \code{newx} must be a matrix, \code{predict} function does not
##'   accept a vector or other formats of \code{newx}.
##' @param s value(s) of the penalty parameter \code{lambda} at which
##'   predictions are required. Default is the entire sequence used to create
##'   the model.
##' @param type type of prediction required. \itemize{ \item Type \code{"link"}
##'   gives the linear predictors for classification problems and gives
##'   predicted response for regression problems. \item Type \code{"class"}
##'   produces the class label corresponding to the maximum probability. Only
##'   available for classification problems.}
##' @param \dots Not used. Other arguments to predict.
##' @return The object returned depends on type.
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'
##' Maintainer: Yi Yang <yi.yang6@mcgill.ca>
##' @seealso \code{\link{coef}} method
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
##' print(predict(m1, type = "class",newx = FHT$x[2:5, ]))
##'
##' @export
predict.gcdnet <- function(object, newx, s = NULL, type = c("class", "link"),
                           ...) NextMethod("predict")

##' Model predictions
##'
##' \code{predict} is a generic function for predictions from the results of
##' various model fitting functions. The function invokes particular
##' \emph{methods} which depend on the \code{\link{class}} of the first
##' argument.
##'
##' @param object a model object for which prediction is desired.
##' @param \dots additional arguments affecting the predictions produced.
##'
##' @return The form of the value returned by \code{predict} depends on the
##'   class of its argument. See the documentation of the particular methods for
##'   details of what is produced by that method.
##'
##' @seealso \code{\link{predict.gcdnet}}, \code{\link{predict.erpath}},
##'   \code{\link{predict.lspath}}, \code{\link{predict.hsvmpath}},
##'   \code{\link{predict.logitpath}}, \code{\link{predict.sqsvmpath}}.
##'
##' @export
predict <- function(object,...) UseMethod("predict")
