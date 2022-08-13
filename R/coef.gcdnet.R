##' Get coefficients or make coefficient predictions from a "gcdnet" object.
##'
##' Computes the coefficients or returns a list of the indices of the nonzero
##' coefficients at the requested values for \code{lambda} from a fitted
##' \code{\link{gcdnet}} object.
##'
##' \code{s} is the new vector at which predictions are requested. If \code{s}
##' is not in the lambda sequence used for fitting the model, the \code{coef}
##' function will use linear interpolation to make predictions. The new values
##' are interpolated using a fraction of coefficients from both left and right
##' \code{lambda} indices.
##'
##' @aliases coef.gcdnet coef.hsvmpath coef.sqsvmpath coef.logitpath coef.lspath
##'   coef.erpath
##' @param object fitted \code{\link{gcdnet}} model object.
##' @param s value(s) of the penalty parameter \code{lambda} at which
##'   predictions are required. Default is the entire sequence used to create
##'   the model.
##' @param type type \code{"coefficients"} computes the coefficients at the
##'   requested values for \code{s}. Type \code{"nonzero"} returns a list of the
##'   indices of the nonzero coefficients for each value of \code{s}. Default is
##'   \code{"coefficients"}.
##' @param \dots not used. Other arguments to predict.
##'
##' @return The object returned depends on type.
##'
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr
##'
##'   Maintainer: Yi Yang <yi.yang6@mcgill.ca>
##'
##' @seealso \code{\link{predict.gcdnet}} method
##'
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
##'
##' @keywords models regression
##'
##' @examples
##'
##' data(FHT)
##' fit1 <- gcdnet(x = FHT$x,y = FHT$y)
##' coef(fit1, type = "coef", s = c(0.1,0.005))
##' coef(fit1, type = "nonzero")
##'
##' @export
coef.gcdnet <- function(object, s = NULL,
                        type = c("coefficients", "nonzero"),
                        ...) NextMethod("coef")

##' Extract Model Coefficients
##'
##' \code{coef} is a generic function which extracts model coefficients from
##' objects returned by modeling functions. \code{coefficients} is an
##' \emph{alias} for it.
##'
##' @param object an object for which the extraction of model coefficients is
##'   meaningful.
##' @param \ldots other arguments.
##'
##' @return Coefficients extracted from the model object \code{object}.
##'
##' @seealso \code{\link{coef.gcdnet}}, \code{\link{coef.erpath}},
##'   \code{\link{coef.lspath}}, \code{\link{coef.hsvmpath}},
##'   \code{\link{coef.logitpath}}, \code{\link{coef.sqsvmpath}}.
##'
##' @export
coef <- function(object, ...) UseMethod("coef")
