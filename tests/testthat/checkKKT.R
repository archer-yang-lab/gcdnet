library(gcdnet)
set.seed(123456)
n=100
p=400
# n=4
# p=8
x=matrix(rnorm(n*p),n,p)
y=rnorm(n)
pf = abs(rnorm(p))
pf2 = abs(rnorm(p))
lambda2=abs(rnorm(1))
eps=1e-12

### include intercept
m1 <-gcdnet(x=x,y=y,eps=eps,pf=pf,pf2=pf2,
            standardize=FALSE,lambda2=lambda2,
            intercept=TRUE, method="ls")

### check KKT
B <- as.matrix(m1$beta)
for (l in 1:length(m1$lambda)) {
  ri=y-(x%*%B[,l]+m1$b0[l])
  for (j in 1:p) {
    L = ri
    xl <- crossprod(x[,j], L)
    if(B[j,l]!=0) {
      AA=-xl/n + lambda2*B[j,l]*pf2[j]+ pf[j]*m1$lambda[l]*sign(B[j,l])
      if(abs(AA)>=sqrt(eps)) print(paste0("lam",l,"-beta",j,":AA=",AA))
    }
    else {
      BB <- abs(-xl/n)-m1$lambda[l]*pf[j]
      if (BB > 0) print(paste("BB=",BB,"-l=",l,"-j=",j,sep=""))
    }
  }
}


### exclude intercept
m2 <-gcdnet(x=x,y=y,eps=eps,pf=pf,pf2=pf2,
            standardize=FALSE,lambda2=lambda2,
            intercept=FALSE, method="ls")

### check KKT
B <- as.matrix(m2$beta)
for (l in 1:length(m2$lambda)) {
  ri=y-(x%*%B[,l]+m2$b0[l])
  for (j in 1:p) {
    L = ri
    xl <- crossprod(x[,j], L)
    if(B[j,l]!=0) {
      AA=-xl/n+lambda2*B[j,l]*pf2[j]+pf[j]*m2$lambda[l]*sign(B[j,l])
      if(abs(AA)>=sqrt(eps)) print(paste0("lam",l,"-beta",j,":AA=",AA))
    }
    else {
      BB <- abs(-xl/n)-m2$lambda[l]*pf[j]
      if (BB > 0) print(paste("BB=",BB,"-l=",l,"-j=",j,sep=""))
    }
  }
}
