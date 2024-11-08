############## General KDE and non-parametric regression #######################
################################################################################

k<-function(x, type="gaussian"){
  if(type=="gaussian"){
    temp<-dnorm(x)
  }
  if(type=="logistic"){
    temp<-(2/pi)*exp(x)/(1+exp(2*x))
  }
  if(type=="epanechnikov"){
    temp<-.75*(1-x^2)
  }
  return(temp)
}

# Non-parametric regression
w<-function(x, h, X, type){
  num<-k((X-x)/h, type)/h
  denom<-sum(k((X-x)/h, type)/h)
  temp<-num/denom
  return(temp)
}

r_x<-function(x, h, X, Y, type){
  weight<-w(x, h, X, type)
  temp<-sum(weight*Y)
  return(temp)
}


r<-function(X, Y, intercept, n=100, h=1/10, type="gaussian"){
  if(is.matrix(X)==FALSE){
    x<-seq(min(X), max(X), length.out=n)
    y_hat<-c()
    for(point in x){
      y_hat<-c(y_hat, r_x(point, h, X, Y, type))
    }
  }
  if(is.matrix(X)==TRUE){
    if(intercept==TRUE){
      d<-ncol(X)-1
      x<-matrix(rep(0,n*d), ncol=d)
      for(j in 1:d){
        x[,j]<-seq(min(X[,j+1]), max(X[,j+1]), length.out=n)
      }
      X<-X[,-1]
    }
    if(intercept==FALSE){
      d<-ncol(X)
      x<-matrix(rep(0,n*d), ncol=d)
      for(j in 1:d){
        x[,j]<-seq(min(X[,j]), max(X[,j]), length.out=n)
      }
    }
    y_hat<-c()
    for(i in 1:n){
      point<-x[i,]
      num<-c()
      for(l in 1:nrow(X)){
        num<-c(num, prod(k((X[l,]-point)/h, type)/h))
      }
      denom<-sum(num)
      W<-num/denom
      y_hat<-c(y_hat,sum(W*Y))
    }
    
  }
  return(list(x=x, y_hat=y_hat))
}


# Normal Example
## Non-parametric regression
n<-100
b<-matrix(c(5.01, 9.22, -12.22, 1, 2), ncol=1)
x1<-seq(1, 2, length.out=n)
x2<-seq(5.01,9.22,length.out=n)
x3<-x1+rgamma(n,1,8)
x4<-x2+rgamma(n,1,4)
X<-matrix(c(rep(1,n),x1^2, x2^2, x3, x4^2), ncol=5)
Y<-X%*%b+rnorm(n, 0, 50)

Y_hat<-r(X, Y, intercept=TRUE, h=c(.75, .5), n=100)
Y_hat$y_hat
sum((Y-Y_hat$y_hat)^2)

plot(x1, Y)
points(sqrt(Y_hat$x[,1]), Y_hat$y_hat, type="l")
plot(x2,Y)
points(sqrt(Y_hat$x[,2]), Y_hat$y_hat, type="l")
plot(x3, Y)
points(Y_hat$x[,3], Y_hat$y_hat, type="l")
plot(x4,Y)
points(sqrt(Y_hat$x[,4]), Y_hat$y_hat, type="l")

# Poisson Example
## Non-parametric regression
n<-100
b<-matrix(c(5.01, 9.22, 12.22), ncol=1)
x1<-seq(1, 2, length.out=n)
x2<-seq(5.01,9.22,length.out=n)
X<-matrix(c(rep(1,n),x1^2, x2), ncol=3)
lam<-exp(X%*%b)
Y<-c()
for(i in 1:n){
  Y<-c(Y, rpois(1, lam[i]))
}

Y_hat<-r(X, Y, intercept=TRUE, h=c(.75, .5), n=100)
plot(x1, Y)
points(sqrt(Y_hat$x[,1]), Y_hat$y_hat, type="l")
plot(x2,Y)
points(Y_hat$x[,2], Y_hat$y_hat, type="l")


