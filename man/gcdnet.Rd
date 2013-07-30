\name{gcdnet}
\alias{gcdnet}
\title{Fits the regularization paths for large margin classifiers}
\description{Fits a regularization path for large margin classifiers at a sequence of regularization parameters lambda.}
\usage{
gcdnet(x, y, nlambda = 100, 
		method = c("hhsvm", "logit", "sqsvm", "ls"),
		lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04), 
		lambda = NULL, lambda2 = 0, 
		pf = rep(1, nvars), pf2 = rep(1, nvars), exclude, 
		dfmax = nvars + 1, pmax = min(dfmax * 1.2, 
	    nvars), standardize = TRUE, eps = 1e-8, maxit = 1e6, delta = 2)
}
\arguments{
		\item{x}{matrix of predictors, of dimension \eqn{N \times p}{N*p}; each row is an observation vector.}

		\item{y}{response variable. This argument should be a two-level factor for classification.}
		\item{nlambda}{the number of \code{lambda} values - default is 100.}
		\item{method}{a character string specifying the loss function to use, valid options are:
		\itemize{
		\item \code{"hhsvm"} Huberized squared hinge loss,
		\item \code{"sqsvm"} Squared hinge loss,
		\item \code{"logit"} logistic loss,
		\item \code{"ls"} least square loss.
		} Default is \code{"hhsvm"}.}
		\item{lambda.factor}{The factor for getting the minimal lambda in \code{lambda} sequence, where \code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}.  \code{max(lambda)} is the smallest value of \code{lambda} for which all coefficients are zero. The default depends on the relationship between \eqn{N} (the number of rows in the matrix of predictors) and \eqn{p}{p} (the number of predictors). If \eqn{N > p}, the default is \code{0.0001},
		close to zero.  If \eqn{N<p}, the default is \code{0.01}.
		A very small value of \code{lambda.factor} will lead to a saturated fit. It takes no effect if there is user-defined \code{lambda} sequence.} 
		\item{lambda}{a user supplied \code{lambda} sequence. Typically, by leaving this option unspecified users can have 
		the program compute its own \code{lambda} sequence based on
		\code{nlambda} and \code{lambda.factor}. Supplying a value of
		\code{lambda} overrides this. It is better to supply
		a decreasing sequence of \code{lambda} values than a single (small) value, if not, the program will sort user-defined \code{lambda} sequence in decreasing order automatically.}

		\item{lambda2}{regularization parameter \eqn{\lambda_2}{lambda2} for the quadratic penalty of the 
		coefficients.}

		\item{pf}{L1 penalty factor of length \eqn{p}{p} used for adaptive LASSO or adaptive elastic net. Separate L1 penalty weights can be applied to each coefficient of \eqn{\beta}{beta} to allow
		differential L1 shrinkage. Can be 0 for some variables, which implies
		no L1 shrinkage, and results in that variable always being included in the
		model. Default is 1 for all variables (and implicitly infinity for
		variables listed in \code{exclude}).}
		
		\item{pf2}{L2 penalty factor of length \eqn{p}{p} used for adaptive LASSO or adaptive elastic net. Separate L2 penalty weights can be applied to each coefficient of \eqn{\beta}{beta} to allow
		differential L2 shrinkage. Can be 0 for some variables, which implies
		no L2 shrinkage. Default is 1 for all variables.}

		\item{exclude}{indices of variables to be excluded from the
		model. Default is none. Equivalent to an infinite penalty factor.}

		\item{dfmax}{limit the maximum number of variables in the
		model. Useful for very large \eqn{p}, if a partial path is desired. Default is \eqn{p+1}.}

		\item{pmax}{limit the maximum number of variables ever to be nonzero. For example once \eqn{\beta} enters the model, no matter how many times it exits or re-enters model through the path, it will be counted only once. Default is \code{min(dfmax*1.2,p)}.}

		\item{standardize}{logical flag for variable standardization, prior to
		fitting the model sequence. If \code{TRUE}, x matrix is normalized such that sum squares of each column \eqn{\sum^N_{i=1}x_{ij}^2/N=1}{<Xj,Xj>/N=1}. Note that x is always centered (i.e. \eqn{\sum^N_{i=1}x_{ij}=0}{sum(Xj)=0}) no matter \code{standardize} is \code{TRUE} or \code{FALSE}. The coefficients are always returned on
		the original scale. Default is is \code{TRUE}.}

		\item{eps}{convergence threshold for coordinate majorization descent. Each inner
		coordinate majorization descent loop continues until the relative change in any
		coefficient (i.e. \eqn{\max_j(\beta_j^{new}-\beta_j^{old})^2}{max(j)(beta_new[j]-beta_old[j])^2}) is less than \code{eps}. For HHSVM i.e. \code{method="hhsvm"}, it is \eqn{\frac{2}{\delta}\max_j(\beta_j^{new}-\beta_j^{old})^2}{2*max(j)(beta_new[j]-beta_old[j])^2/delta}. Defaults value is \code{1e-8}.}

		\item{maxit}{maximum number of outer-loop iterations allowed at fixed lambda value. Default is 1e6. If models do not converge, consider increasing \code{maxit}.}

		\item{delta}{the parameter \eqn{\delta}{delta} in HHSVM model. Default is 2.}
}

\details{
Note that the objective function in \code{gcdnet} is
\deqn{Loss(y, X, beta))/N + lambda1 * |beta| + 0.5 * lambda2 * beta^2}
where the penalty is a combination of L1 and L2 term. Users can specify the loss function to use, options include Huberized squared hinge loss, Squared hinge loss, least square loss and logistic regression. Users can also tweak the penalty by choosing different \eqn{lambda2} and penalty factor. 

For computing speed reason, if models are not converging or running slow, consider increasing \code{eps}, decreasing
\code{nlambda}, or increasing \code{lambda.factor} before increasing
\code{maxit}.

\strong{FAQ:}

\bold{Question: }\dQuote{\emph{I couldn't get an idea how to specify an option to get adaptive LASSO, how to specify an option to get elastic net and adaptive elastic net? Could you please give me a quick hint?}}

\bold{Answer: } \code{lambda2} is the regularize parameter for L2 penalty part. To use LASSO, set \code{lambda2=0}. To use elastic net, set \code{lambda2} as nonzero.

\code{pf} is the L1 penalty factor of length \eqn{p}{p} (\eqn{p}{p} is the number of predictors). Separate L1 penalty weights can be applied to each coefficient to allow differential L1 shrinkage. Similiarly \code{pf2} is the L2 penalty factor of length \eqn{p}{p}.

To use adaptive LASSO, you should set \code{lambda2=0} and also specify \code{pf} and \code{pf2}. To use adaptive elastic net, you should set \code{lambda2} as nonzero and specify \code{pf} and \code{pf2},  

For example

\preformatted{
library('gcdnet')

# Dataset N = 100, p = 10
x_log <- matrix(rnorm(100*10),100,10)
y_log <- sample(c(-1,1),100,replace=TRUE)

# LASSO
m <- gcdnet(x=x_log,y=y_log,lambda2=0,method="log")
plot(m)

# elastic net with lambda2 = 1 
m <- gcdnet(x=x_log,y=y_log,lambda2=1,method="log")
plot(m)

# adaptive lasso with penalty factor 
# pf = 0.5 0.5 0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0
m <- gcdnet(x=x_log,y=y_log,lambda2=0,method="log",
pf=c(rep(0.5,5),rep(1,5)))
plot(m)

# adaptive elastic net with lambda2 = 1 and penalty factor pf = c(rep(0.5,5),rep(1,5))
# pf2 = 3 3 3 3 3 1 1 1 1 1

m <- gcdnet(x=x_log,y=y_log,lambda2=1,method="log",
pf=c(rep(0.5,5),rep(1,5)), 
pf2 = c(rep(3,5),rep(1,5)))
plot(m)
}

\bold{Question: }\dQuote{\emph{what is the meaning of the parameter \code{pf}? On the package documentation, it said \code{pf} is the penalty weight applied to each coefficient of beta?}}

\bold{Answer: } Yes, \code{pf} and \code{pf2} are L1 and L2 penalty factor of length \eqn{p}{p} used for adaptive LASSO or adaptive elastic net. 0 means that the feature (variable) is always excluded, 1 means that the feature (variable) is included with weight 1.

\bold{Question: }\dQuote{\emph{Why does predict function not work? predict should return the predicted probability of the positive class. Instead I get:}}
\preformatted{
Error in as.matrix(as.matrix(cbind2(1, newx)) %*% nbeta) : 
  error in evaluating the argument 'x' in selecting 
a method for function 'as.matrix': Error in t(.Call(Csparse_dense_crossprod, y,
 t(x))) : 
  error in evaluating the argument 'x' in selecting 
a method for function 't': Error: Cholmod error 'X and/or Y have wrong dimensions' 
at file ../MatrixOps/cholmod_sdmult.c, line 90?
}

\dQuote{\emph{Using the Arcene dataset and executing the following code will give the above error:}}
\preformatted{
library(gcdnet)
arc <- read.csv("arcene.csv", header=FALSE)
fit <- gcdnet(arc[,-10001], arc[,10001], standardize=FALSE, method="logit")
pred <- rnorm(10000)
predict(fit, pred, type="link")
}

\bold{Answer: } It is actually NOT a bug of gcdnet. When make prediction using a new matrix x, each observation of x should be arranged as a row of a matrix. In your code, because "pred" is a vector, you need to convert "pred" into a matrix, try the following code:
\preformatted{
pred <- rnorm(10000)
pred <- matrix(pred,1,10000)
predict(fit, pred, type="link")
}
}


\value{
An object with S3 class \code{\link{gcdnet}}.
		\item{call}{the call that produced this object}
		\item{b0}{intercept sequence of length \code{length(lambda)}}
		\item{beta}{a \code{p*length(lambda)} matrix of coefficients, stored as a sparse matrix (\code{dgCMatrix} class, the standard class for sparse numeric matrices in the \code{Matrix} package.). To convert it into normal type matrix use \code{as.matrix()}.}
		\item{lambda}{the actual sequence of \code{lambda} values used}
		\item{df}{the number of nonzero coefficients for each value of
		\code{lambda}.}
		\item{dim}{dimension of coefficient matrix (ices)}
		\item{npasses}{total number of iterations (the most inner loop) summed over all lambda values}
		\item{jerr}{error flag, for warnings and errors, 0 if no error.}
}

\author{Yi Yang and Hui Zou\cr
Maintainer: Yi Yang  <yiyang@umn.edu>}
\references{
Yang, Y. and Zou, H. (2012), "An Efficient Algorithm for Computing The HHSVM and Its Generalizations," \emph{Journal of Computational and Graphical Statistics}, 22, 396-415.\cr
BugReport: \url{http://code.google.com/p/gcdnet/}\cr
}


\seealso{\code{plot.gcdnet}}
\examples{
data(FHT)
# solution paths for the least squares
m0 <- gcdnet(x=FHT$x,y=FHT$y_reg,lambda2=1,method="ls")
plot(m0)
# solution paths for the HHSVM
m1 <- gcdnet(x=FHT$x,y=FHT$y,delta=1,lambda2=1,method="hhsvm")
plot(m1)
# solution paths for the penalized SVM with the squared hinge loss
m2 <- gcdnet(x=FHT$x,y=FHT$y,lambda2=0.1,method="sqsvm")
plot(m2)
# solution paths for the penalized logistic regression
m3 <- gcdnet(x=FHT$x,y=FHT$y,lambda2=0.01,method="logit")
plot(m3)
}
\keyword{models}
\keyword{regression}
