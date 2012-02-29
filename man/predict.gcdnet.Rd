\name{predict.gcdnet}
\alias{predict.gcdnet}
\alias{predict.hsvmpath}
\alias{predict.sqsvmpath}
\alias{predict.logitpath}
\title{make predictions from a "gcdnet" object.}

\description{
Similar to other predict methods, this functions predicts fitted values and class labels from a fitted \code{\link{gcdnet}} object.
}

\usage{
\method{predict}{gcdnet}(object, newx, s = NULL,
type=c("class","link"), ...)
}

\arguments{
		\item{object}{fitted \code{\link{gcdnet}} model object.}
		\item{newx}{matrix of new values for \code{x} at which predictions are
		to be made. Must be a matrix.}
		\item{s}{value(s) of the penalty parameter \code{lambda} at which
		predictions are required. Default is the entire sequence used to
		create the model. }
		\item{type}{
		type of prediction required. 
		\itemize{
				\item Type \code{"link"} gives the
				linear predictors for classification problems.
				\item Type \code{"class"} produces the class label corresponding to the maximum probability.}
		}
		\item{\dots}{Not used. Other arguments to predict. } }
\details{\code{s} is the new vector at which predictions are requested. If \code{s} is not in the lambda sequence used for fitting the model, the \code{predict} function will use linear interpolation to make predictions. The new values are interpolated using a fraction of predicted values from both left and right \code{lambda} indices. }
\value{The object returned depends on type.}


\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yi@stat.umn.edu>}
\references{
 Yi Yang and Hui Zou. \dQuote{An Efficient Algorithm for Computing HHSVM and Its Generalization}.   School of Statistics, University of Minnesota.

\url{http://code.google.com/p/gcdnet/}}

\seealso{\code{\link{coef}} method}
\examples{
data(FHT)
m1 = gcdnet(x=FHT$x,y=FHT$y)
print(predict(m1,type="class",newx=FHT$x[2:5,]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
