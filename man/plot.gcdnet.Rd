\name{plot.gcdnet}
\alias{plot.gcdnet}
\title{Plot coefficients from a "gcdnet" object}
\description{
Produces a coefficient profile plot of the coefficient paths for a
fitted \code{\link{gcdnet}} object.
}
\usage{
\method{plot}{gcdnet}(x, xvar = c("norm", "lambda"), color = FALSE, label = FALSE, ...)
}
\arguments{
		\item{x}{fitted \code{\link{gcdnet}} model}
		\item{xvar}{what is on the X-axis. \code{"norm"} plots against the
		L1-norm of the coefficients, \code{"lambda"} against the log-lambda
		sequence.}
		\item{color}{if \code{TRUE}, plot the curves with rainbow colors. \code{FALSE} is gray colors. Default is \code{FALSE}}
		\item{label}{if \code{TRUE}, label the curves with variable sequence numbers. Default is \code{FALSE}}
		\item{\dots}{other graphical parameters to plot}
}
\details{
A coefficient profile plot is produced.
}
\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yi@stat.umn.edu>}
\references{
 Yi Yang and Hui Zou. \dQuote{An Efficient Algorithm for Computing HHSVM and Its Generalization}.   School of Statistics, University of Minnesota.

\url{http://code.google.com/p/gcdnet/}}

\examples{
data(FHT)
m1 = gcdnet(x=FHT$x,y=FHT$y)
par(mfrow=c(1,3))
plot(m1) # plots against the L1-norm of the coefficients
plot(m1,xvar="lambda",label=TRUE) # plots against the log-lambda sequence
plot(m1,color=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
