


---


# R package: gcdnet #



Author: Yi Yang and Hui Zou



Reference: Yang, Y. and Zou, H. (2012), _An Efficient Algorithm for Computing HHSVM and Its Generalization_, Journal of Computational and Graphical Statistics. Accepted.



This package is Copyright (C) 2012 Yi Yang and Hui Zou.  All rights reserved by the authors.



---


# Overview #



This an [R](http://cran.r-project.org) package implements a generalized coordinate descent (GCD) algorithm for computing the solution path of the hybrid Huberized support vector machine (HHSVM) and its generalization, including the elastic net penalized SVM with the squared hinge loss and the elastic net penalized logistic regression.



---

# Installation #




## Install from CRAN ##

The package is now available on the [CRAN](http://cran.r-project.org) gcdnet site: http://cran.r-project.org/web/packages/gcdnet.

It can be installed using the commands
```
install.packages("gcdnet")
library(gcdnet)
```

## Install from this Google Code site ##

### R under Windows ###

  1. Go to the "Downloads" tab on this site and download the most recent zip-file (Windows binary) to your hard drive.
  1. Open R and select from the "Packages" menu at the top of the R console "Install package(s) from local zip files...".
  1. Navigate to the folder where you saved the zip-file, select the zip-file and click open.

### R under Linux / Mac OS X ###

  1. Go to the "Downloads" tab on this site and download the most recent tar.gz-file (source code) to your hard drive.
  1. Open a terminal window and navigate to the folder where you saved the tar.gz-file and type: R CMD INSTALL
  1. gcdnet\_x.x.x.tar.gz (substitute x.x.x with the version number)



---


# Documentation #




Download the full pdf reference [documentation](http://code.google.com/p/gcdnet/downloads/list) which contains a function index



---


# Getting Started #




After you installed the package type into R:
```
library(gcdnet)
```
Type
```
?gcdnet
```
to see overall documentation.

See some example here

```

# 0. Load gcdnet library and FHT data.
library(gcdnet)
data(FHT)

# 1. fit the solution paths for the HHSVM
m1 <- gcdnet(x=FHT$x,y=FHT$y,delta=1,lambda2=1, method="hhsvm")

# 2. make a plot for solution paths.
plot(m1)

# 3. fit the solution paths for the penalized SVM with the squared hinge loss.
m2 <- gcdnet(x=FHT$x, y=FHT$y, lambda2=0.1, method="sqsvm")

# 4. fit the solution paths for the penalized logistic regression.
m3 <- gcdnet(x=FHT$x, y=FHT$y, lambda2=0.01, method="logit")

# 5. Does 5-fold cross-validation for gcdnet and returns a optimal value for lambda (default model is "hhsvm")
cv=cv.gcdnet(FHT$x, FHT$y, nfolds=5, lambda2=0.1)

# 6. plot the cross-validation curve produced by 'cv.gcdnet()'
plot(cv)
```


---


# FAQ #

## Q: After quickly going through the manual, I couldn't get an idea how to specify an option to get adaptive LASSO, how to specify an option to get elastic net and adaptive elastic net? Could you please give me a quick hint? ##

## A: ##
"lambda2" is the regularize parameter for L2 penalty part. To use LASSO, set "lambda2=0" To use elastic net, set "lambda2" as nonzero.

"pf" is the L1 penalty factor of length p. Separate L1 penalty weights can be applied to each coefficient of β to allow differential L1 shrinkage. Can be 0 for some variables, which implies no L1 shrinkage, and results in that variable always being included in the model. Default is 1 for all variables .

Similiarly "pf2" is the L2 penalty factor of length p.

To use adaptive LASSO, you can specify weight to L1 penalty factor "pf", and set lambda2=0
To use adaptive elastic net, you can specify both pf and pf2, and set lambda2 as nonzero.

For example

```
library('gcdnet')

# Dataset N = 100, p = 10
x_ls <- matrix(rnorm(100*10),100,10)
y_ls <- rnorm(100)

# LASSO
m0 <- gcdnet(x=x_ls,y=y_ls,lambda2=0,method="ls")
plot(m0)

# elastic net with lambda2 = 1 
m0 <- gcdnet(x=x_ls,y=y_ls,lambda2=1,method="ls")
plot(m0)

# adaptive lasso with penalty factor pf = 0.5 0.5 0.5 0.5 0.5 1.0 1.0 1.0 1.0 1.0
m0 <- gcdnet(x=x_ls,y=y_ls,lambda2=0,method="ls",pf=c(rep(0.5,5),rep(1,5)))
plot(m0)

# adaptive elastic net with lambda2 = 1 and penalty factor pf = c(rep(0.5,5),rep(1,5))
# pf2 = 3 3 3 3 3 1 1 1 1 1

m0 <- gcdnet(x=x_ls,y=y_ls,lambda2=1,method="ls",pf=c(rep(0.5,5),rep(1,5)), pf2 = c(rep(3,5),rep(1,5)))
plot(m0)
```

## Q: what is the meaning of the parameter "pf"? On the package documentation, it said "pf" is the penalty weight applied to each coefficient of beta? ##

## A: ##
Yes, "pf" and "pf2" are L1 and L2 penalty factor of length p used for adaptive LASSO or adaptive elastic net. 0 means that the feature (variable) is always excluded, 1 means that the feature (variable) is included with weight 1.

# Screen Shot #

![http://gcdnet.googlecode.com/files/Screenshot1.jpg](http://gcdnet.googlecode.com/files/Screenshot1.jpg)
![http://gcdnet.googlecode.com/files/Screenshot2.jpg](http://gcdnet.googlecode.com/files/Screenshot2.jpg)
![http://gcdnet.googlecode.com/files/Screenshot3.jpg](http://gcdnet.googlecode.com/files/Screenshot3.jpg)