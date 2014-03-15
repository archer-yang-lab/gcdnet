library(gcdnet)
data(FHT)
m0 <- gcdnet(x=FHT$x,y=FHT$y_reg,lambda2=1,method="ls")#
plot(m0)#
#
# 2. solution paths for the HHSVM#
m1 <- gcdnet(x=FHT$x,y=FHT$y,delta=1,lambda2=1,method="hhsvm")#
plot(m1)#
#
# 3. solution paths for the penalized SVM with the squared hinge loss#
m2 <- gcdnet(x=FHT$x,y=FHT$y,lambda2=0.1,method="sqsvm")#
plot(m2)
# 4. solution paths for the adaptive elastic net penalized logistic regression#
#
p <- ncol(FHT$x)#
# set the first three L1 penalty weights as 2 and the rest are 1#
pf = c(2,2,2,rep(1,p-3))#
#
# set the last three L2 penalty weights as 0.5 and the rest are 0.5#
pf2 = c(rep(1,p-3),0.5,0.5,0.5)
pf2 = c(rep(1,p-3),0.5,0.5,0.5)#
m3 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,pf2=pf2,lambda2=0.01,method="logit")#
plot(m3)
# 1. solution paths for the LASSO penalized least squares#
# to use LASSO set lambda2 = 0#
#
m1 <- gcdnet(x=FHT$x,y=FHT$y_reg,lambda2=0,method="ls")#
plot(m1)
# 2. solution paths for the HHSVM#
#
m2 <- gcdnet(x=FHT$x,y=FHT$y,delta=1,lambda2=1,method="hhsvm")#
plot(m2)
# set the first three L1 penalty weights as 10 and the rest are 1#
pf = c(10,10,10,rep(1,p-3))#
m3 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,lambda2=0.1,method="sqsvm")#
plot(m3)
# set the first three L1 penalty weights as 10 and the rest are 1#
pf = c(1,1,1,rep(0.01,p-3))
m3 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,lambda2=0.1,method="sqsvm")
plot(m3)
p <- ncol(FHT$x)#
# set the first three L1 penalty weights as 10 and the rest are 1#
pf = c(1,1,1,rep(0.01,p-3))#
m3 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,lambda2=0.1,method="sqsvm")#
plot(m3)
p <- ncol(FHT$x)#
# set the first three L1 penalty weights as 10 and the rest are 1#
pf = c(0.01,0.01,0.01,rep(1,p-3))#
m3 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,lambda2=0.1,method="sqsvm")#
plot(m3)
pf = c(0.1,0.1,0.1,rep(1,p-3))#
m3 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,lambda2=0.1,method="sqsvm")#
plot(m3)
p <- ncol(FHT$x)#
# set the first three L1 penalty weights as 10 and the rest are 1#
pf = c(0.1,0.1,0.1,rep(1,p-3))#
m3 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,lambda2=0.1,method="sqsvm")#
plot(m3)
setwd('/Users/emeryyi/Dropbox/Research/googleproject/gcdnet/man')
# 3. solution paths for the adaptive LASSO penalized SVM #
# with the squared hinge loss. To use the adaptive LASSO, #
# set lambda2 = 0 and meanwhile specify the L1 penalty weights.#
p <- ncol(FHT$x)#
# set the first three L1 penalty weights as 0.1 and the rest are 1#
pf = c(0.1,0.1,0.1,rep(1,p-3))#
m3 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,lambda2=0.1,method="sqsvm")#
plot(m3)#
#
# 4. solution paths for the adaptive elastic net penalized logistic regression#
#
p <- ncol(FHT$x)#
# set the first three L1 penalty weights as 10 and the rest are 1#
pf = c(10,10,10,rep(1,p-3))#
# set the last three L2 penalty weights as 0.1 and the rest are 1#
pf2 = c(rep(1,p-3),0.1,0.1,0.1)#
# set the L2 penalty parameter lambda2=0.01#
m4 <- gcdnet(x=FHT$x,y=FHT$y,pf=pf,pf2=pf2,lambda2=0.01,method="logit")#
plot(m4)
