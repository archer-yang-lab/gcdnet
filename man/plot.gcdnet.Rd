% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.gcdnet.R
\name{plot.gcdnet}
\alias{plot.gcdnet}
\title{Plot coefficients from a "gcdnet" object}
\usage{
\method{plot}{gcdnet}(x, xvar = c("norm", "lambda"), color = FALSE, label = FALSE, ...)
}
\arguments{
\item{x}{fitted \code{\link{gcdnet}} model}

\item{xvar}{what is on the X-axis. \code{"norm"} plots against the L1-norm
of the coefficients, \code{"lambda"} against the log-lambda sequence.}

\item{color}{if \code{TRUE}, plot the curves with rainbow colors.
\code{FALSE} is gray colors. Default is \code{FALSE}}

\item{label}{if \code{TRUE}, label the curves with variable sequence
numbers. Default is \code{FALSE}}

\item{\dots}{other graphical parameters to plot}
}
\description{
Produces a coefficient profile plot of the coefficient paths for a fitted
\code{\link{gcdnet}} object. This function is modified based on the
\code{plot} function from the \code{glmnet} package.
}
\details{
A coefficient profile plot is produced.
}
\examples{

data(FHT)
m1 <- gcdnet(x = FHT$x,y = FHT$y)
par(mfrow = c(1,3))
plot(m1) # plots against the L1-norm of the coefficients
plot(m1,xvar = "lambda",label = TRUE) # plots against the log-lambda sequence
plot(m1,color = TRUE)

}
\references{
Yang, Y. and Zou, H. (2012).
  "An Efficient Algorithm for Computing The HHSVM and Its Generalizations."
  \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
  BugReport: \url{https://github.com/emeryyi/gcdnet}\cr

  Gu, Y., and Zou, H. (2016).
  "High-dimensional generalizations of asymmetric least squares regression and their applications."
  \emph{The Annals of Statistics}, 44(6), 2661–2694.\cr

  Friedman, J., Hastie, T., and Tibshirani, R. (2010).
  "Regularization paths for generalized linear models via coordinate descent."
  \emph{Journal of Statistical Software, 33, 1.}\cr
  \url{https://www.jstatsoft.org/v33/i01/}
}
\author{
Yi Yang, Yuwen Gu and Hui Zou\cr Maintainer: Yi Yang
  <yi.yang6@mcgill.ca>
}
\keyword{models}
\keyword{regression}
