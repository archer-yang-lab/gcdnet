\name{plot.cv.gcdnet}
\alias{plot.cv.gcdnet}
\title{plot the cross-validation curve produced by cv.gcdnet}
\description{
Plots the cross-validation curve, and upper and lower standard deviation
curves, as a function of the \code{lambda} values used.}
\usage{
\method{plot}{cv.gcdnet}(x, sign.lambda, ...)
}
\arguments{
	\item{x}{fitted \code{\link{cv.gcdnet}} object}
	\item{sign.lambda}{either plot against \code{log(lambda)} (default) or
	its negative if \code{sign.lambda=-1}.}
	\item{\dots}{other graphical parameters to plot}
}
\details{
A plot is produced.}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yi@stat.umn.edu>}
\references{
 Yi Yang and Hui Zou. \dQuote{An Efficient Algorithm for Computing HHSVM and Its Generalization}.   School of Statistics, University of Minnesota.

\url{http://code.google.com/p/gcdnet/}}

\seealso{\code{\link{cv.gcdnet}}.}
\examples{
data(colon)
set.seed(2011)
cv=cv.gcdnet(colon$x, colon$y,
lambda2 = 1, pred.loss="loss", nfolds=5)
plot(cv)
}
\keyword{models}
\keyword{regression}
