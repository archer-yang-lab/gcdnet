\name{cv.gcdnet}
\alias{cv.gcdnet}
\alias{cv.hsvmpath}
\alias{cv.sqsvmpath}
\alias{cv.logitpath}
\alias{cv.lspath}
\alias{cv.erpath}
\title{Cross-validation for gcdnet}
\description{Does k-fold cross-validation for gcdnet, produces a plot,
and returns a value for \code{lambda}. This function is modified based on the \code{cv} function from the \code{glmnet} package.}
\usage{
cv.gcdnet(x, y, lambda, pred.loss, nfolds, foldid, delta, omega,...)
}
\arguments{
		\item{x}{\code{x} matrix as in \code{\link{gcdnet}}.}
		\item{y}{response variable or class label \code{y} as in \code{\link{gcdnet}}.}
		\item{lambda}{optional user-supplied lambda sequence; default is
		\code{NULL}, and \code{\link{gcdnet}} chooses its own sequence.}
		\item{nfolds}{number of folds - default is 5. Although \code{nfolds}
		can be as large as the sample size (leave-one-out CV), it is not
		recommended for large datasets. Smallest value allowable is \code{nfolds=3}.}
		\item{foldid}{an optional vector of values between 1 and \code{nfold}
		identifying what fold each observation is in. If supplied,
		\code{nfold} can be missing.}
		\item{pred.loss}{loss function to use for cross-validation error. Valid options are:
		\itemize{
		\item \code{"loss"} Margin based loss function. When use least square loss \code{"ls"}, it gives mean square error (MSE). When use expectile regression loss \code{"er"}, it gives asymmetric mean square error (AMSE).
		\item \code{"misclass"} only available for classification: it gives misclassification error. }
		Default is \code{"loss"}.}
		\item{delta}{parameter \eqn{\delta}{delta} only used in HHSVM for computing margin based loss function, only available for \code{pred.loss = "loss"}.}
		\item{omega}{parameter \eqn{\omega}{omega} only used in expectile regression. Only available for \code{pred.loss = "loss"}.}
		\item{\dots}{other arguments that can be passed to gcdnet.}

}
\details{The function runs \code{\link{gcdnet}} \code{nfolds}+1 times; the
first to get the \code{lambda} sequence, and then the remainder to
compute the fit with each of the folds omitted. The average error and standard deviation over the
folds are computed. 
}
\value{an object of class \code{\link{cv.gcdnet}} is returned, which is a
list with the ingredients of the cross-validation fit.
		\item{lambda}{the values of \code{lambda} used in the fits.}
		\item{cvm}{the mean cross-validated error - a vector of length
		\code{length(lambda)}.}
		\item{cvsd}{estimate of standard error of \code{cvm}.}
		\item{cvupper}{upper curve = \code{cvm+cvsd}.}
		\item{cvlower}{lower curve = \code{cvm-cvsd}.}
		\item{nzero}{number of non-zero coefficients at each \code{lambda}.}
		\item{name}{a text string indicating type of measure (for plotting
		purposes).}
		\item{gcdnet.fit}{a fitted \code{\link{gcdnet}} object for the full data.}
		\item{lambda.min}{The optimal value of \code{lambda} that gives minimum
		cross validation error \code{cvm}.}
		\item{lambda.1se}{The largest value of \code{lambda} such that error is
		within 1 standard error of the minimum.}
}

\author{Yi Yang, Yuwen Gu and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{https://github.com/emeryyi/fastcox.git}\cr


Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}\cr
\url{http://www.jstatsoft.org/v33/i01/}}


\seealso{\code{\link{gcdnet}}, \code{\link{plot.cv.gcdnet}}, \code{\link{predict.cv.gcdnet}}, and \code{\link{coef.cv.gcdnet}} methods.}
\examples{
# fit an elastic net penalized HHSVM 
# with lambda2 = 0.1 for the L2 penalty. Use the 
# misclassification rate as the cross validation 
# prediction loss. Use five-fold CV to choose 
# the optimal lambda for the L1 penalty.

data(FHT)
set.seed(2011)
cv=cv.gcdnet(FHT$x, FHT$y, method ="hhsvm",
lambda2=0.1, pred.loss="misclass", nfolds=5, delta=1.5)
plot(cv)

# fit an elastic net penalized least squares 
# with lambda2 = 0.1 for the L2 penalty. Use the 
# least square loss as the cross validation 
# prediction loss. Use five-fold CV to choose 
# the optimal lambda for the L1 penalty.

set.seed(2011)
cv1=cv.gcdnet(FHT$x, FHT$y_reg, method ="ls",
lambda2=0.1,pred.loss="loss", nfolds=5)
plot(cv1)

# To fit a LASSO penalized logistic regression
# we set lambda2 = 0 to disable the L2 penalty. Use the 
# logistic loss as the cross validation 
# prediction loss. Use five-fold CV to choose 
# the optimal lambda for the L1 penalty.

set.seed(2011)
cv2=cv.gcdnet(FHT$x, FHT$y, method ="logit",
lambda2 = 0, pred.loss="loss", nfolds=5)
plot(cv2)
}
\keyword{models}
\keyword{regression}

