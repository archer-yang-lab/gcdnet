##' Cross-validation for gcdnet
##'
##' Does k-fold cross-validation for gcdnet, produces a plot, and returns a
##' value for \code{lambda}. This function is modified based on the \code{cv}
##' function from the \code{glmnet} package.
##'
##' The function runs \code{\link{gcdnet}} \code{nfolds}+1 times; the first to
##' get the \code{lambda} sequence, and then the remainder to compute the fit
##' with each of the folds omitted. The average error and standard deviation
##' over the folds are computed.
##'
##' @aliases cv.gcdnet cv.hsvmpath cv.sqsvmpath cv.logitpath cv.lspath cv.erpath
##'
##' @param x \code{x} matrix as in \code{\link{gcdnet}}.
##' @param y response variable or class label \code{y} as in
##'   \code{\link{gcdnet}}.
##' @param lambda optional user-supplied lambda sequence; default is
##'   \code{NULL}, and \code{\link{gcdnet}} chooses its own sequence.
##' @param nfolds number of folds - default is 5. Although \code{nfolds} can be
##'   as large as the sample size (leave-one-out CV), it is not recommended for
##'   large datasets. Smallest value allowable is \code{nfolds=3}.
##' @param foldid an optional vector of values between 1 and \code{nfold}
##'   identifying what fold each observation is in. If supplied, \code{nfold}
##'   can be missing.
##' @param pred.loss loss function to use for cross-validation error. Valid
##'   options are: \itemize{ \item \code{"loss"} Margin based loss function.
##'   When use least square loss \code{"ls"}, it gives mean square error (MSE).
##'   When use expectile regression loss \code{"er"}, it gives asymmetric mean
##'   square error (AMSE). \item \code{"misclass"} only available for
##'   classification: it gives misclassification error. } Default is
##'   \code{"loss"}.
##' @param delta parameter \eqn{\delta}{delta} only used in HHSVM for computing
##'   margin based loss function, only available for \code{pred.loss = "loss"}.
##' @param omega parameter \eqn{\omega}{omega} only used in expectile
##'   regression. Only available for \code{pred.loss = "loss"}.
##' @param \dots other arguments that can be passed to gcdnet.
##' @return an object of class \code{\link{cv.gcdnet}} is returned, which is a
##'   list with the ingredients of the cross-validation fit. \item{lambda}{the
##'   values of \code{lambda} used in the fits.} \item{cvm}{the mean
##'   cross-validated error - a vector of length \code{length(lambda)}.}
##'   \item{cvsd}{estimate of standard error of \code{cvm}.}
##'   \item{cvupper}{upper curve = \code{cvm+cvsd}.} \item{cvlower}{lower curve
##'   = \code{cvm-cvsd}.} \item{nzero}{number of non-zero coefficients at each
##'   \code{lambda}.} \item{name}{a text string indicating type of measure (for
##'   plotting purposes).} \item{gcdnet.fit}{a fitted \code{\link{gcdnet}}
##'   object for the full data.} \item{lambda.min}{The optimal value of
##'   \code{lambda} that gives minimum cross validation error \code{cvm}.}
##'   \item{lambda.1se}{The largest value of \code{lambda} such that error is
##'   within 1 standard error of the minimum.}
##' @author Yi Yang, Yuwen Gu and Hui Zou\cr Maintainer: Yi Yang
##'   <yi.yang6@mcgill.ca>
##' @seealso \code{\link{gcdnet}}, \code{\link{plot.cv.gcdnet}},
##'   \code{\link{predict.cv.gcdnet}}, and \code{\link{coef.cv.gcdnet}} methods.
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
##' # fit an elastic net penalized HHSVM with lambda2 = 0.1 for the L2 penalty.
##' # Use the misclassification rate as the cross validation prediction loss.
##' # Use five-fold CV to choose the optimal lambda for the L1 penalty.
##'
##' data(FHT)
##' set.seed(2011)
##' cv <- cv.gcdnet(FHT$x, FHT$y, method = "hhsvm",
##'                 lambda2 = 0.1, pred.loss = "misclass",
##'                 nfolds = 5, delta = 1.5)
##' plot(cv)
##'
##' # fit an elastic net penalized least squares
##' # with lambda2 = 0.1 for the L2 penalty. Use the
##' # least square loss as the cross validation
##' # prediction loss. Use five-fold CV to choose
##' # the optimal lambda for the L1 penalty.
##'
##' set.seed(2011)
##' cv1 <- cv.gcdnet(FHT$x, FHT$y_reg, method ="ls",
##'                  lambda2 = 0.1, pred.loss = "loss",
##'                  nfolds = 5)
##' plot(cv1)
##'
##' # To fit a LASSO penalized logistic regression
##' # we set lambda2 = 0 to disable the L2 penalty. Use the
##' # logistic loss as the cross validation
##' # prediction loss. Use five-fold CV to choose
##' # the optimal lambda for the L1 penalty.
##'
##' set.seed(2011)
##' cv2 <- cv.gcdnet(FHT$x, FHT$y, method ="logit",
##'                  lambda2 = 0, pred.loss="loss",
##'                  nfolds=5)
##' plot(cv2)
##'
##' @export
cv.gcdnet <- function(x, y, lambda = NULL, pred.loss = c("misclass", "loss"),
                      nfolds = 5, foldid, delta = 2, omega = 0.5, ...) {
  if (missing(pred.loss))
    pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
  N <- nrow(x)
  ## Fit the model once to get dimensions etc of output
  y <- drop(y)
  gcdnet.object <- gcdnet(x, y, lambda = lambda, delta = delta,
                          omega = omega, ...)
  lambda <- gcdnet.object$lambda
  ## predict -> coef
  nz <- sapply(coef(gcdnet.object, type = "nonzero"), length)
  if (missing(foldid))
    foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
  if (nfolds < 3)
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  outlist <- as.list(seq(nfolds))
  ## Now fit the nfold models and store them
  for (i in seq(nfolds)) {
    which <- foldid == i
    y_sub <- y[!which]
    outlist[[i]] <- gcdnet(x = x[!which, , drop = FALSE],
                           y = y_sub, lambda = lambda,
                           delta = delta, omega = omega, ...)
  }
  ## What to do depends on the pred.loss and the model fit
  fun <- paste("cv", class(gcdnet.object)[[2]], sep = ".")
  cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid,
                               pred.loss, delta, omega))
  cvm <- cvstuff$cvm
  cvsd <- cvstuff$cvsd
  cvname <- cvstuff$name
  out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd,
              cvupper = cvm + cvsd, cvlo = cvm - cvsd,
              nzero = nz, name = cvname, gcdnet.fit = gcdnet.object)
  lamin <- getmin(lambda, cvm, cvsd)
  obj <- c(out, as.list(lamin))
  class(obj) <- "cv.gcdnet"
  obj
}
