\name{FHT}
\alias{FHT}
\docType{data}
\title{FHT data introduced in Friedman et al. (2010).}
\description{
The \code{FHT} data set has n = 50 observations and p = 100 predictors. The covariance between predictors Xj and Xj' has the same correlation 0.5. See details in Friedman et al. (2010).
}
\usage{data(FHT)}
\format{
  This data frame contains the following columns:
  \describe{
    \item{x1}{a matrix with 100 rows and 5000 columns}
    \item{y1}{class labels}
  }
}
\references{
Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths for generalized
linear models via coordinate descent," \emph{Journal of Statistical Software, 33, 1.}
}
\examples{
data(FHT)
}
\keyword{datasets}
