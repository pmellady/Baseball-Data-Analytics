############### Non-Parametric Statistics: Kernel Based Methods ################
### This code examines kernel based regression and kernel density estimation ###
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


r<-function(X, Y, n=100, h=1/10, type="gaussian"){
  n<-length(X)
  x<-seq(min(X), max(X), length.out=n)
  y_hat<-c()
  for(point in x){
    y_hat<-c(y_hat, r_x(point, h, X, Y, type))
  }
  return(list(x=x, y_hat=y_hat))
}

# Density estimation
p_x<-function(x, h, X, type){
  n<-length(X)
  temp<-sum(k((x-X)/h, type)/h)/n
}

p<-function(X, n=1000, h=1/10, type="gaussian"){
  x<-seq(-10*max(abs(X)),10*max(abs(X)), length.out=n)
  p_hat<-c()
  for(point in x){
    p_hat<-c(p_hat, p_x(point, h, X, type))
  }
  return(list(x=x, p_hat=p_hat))
}

# Normal Example
## Non-parametric regression
n<-100
X<-seq(1, 20, length.out=n)
Y<-sqrt(X)+X*exp(-X)*sin(X)+rnorm(n, 0, .1)


Y_hat<-r(X, Y, h=1/2)
plot(X, Y)
points(Y_hat$x, Y_hat$y_hat, type="l")

## KDE
p_hat<-p(Y, h=1/2)
plot(density(Y), col="orange")
points(p_hat$x, p_hat$p_hat, type="l", col="blue")


# Logistic Example
## Non-parametric regression
n<-100
h<-1/10
b<-5.01
x<-seq(-2, 2, length.out=n)
pi<-exp(x*b)/(1+exp(x*b))
Y<-c()
for(i in 1:n){
  Y<-c(Y, rbinom(1,1,pi[i]))
}

Y_hat<-r(x, Y, h=1/2)
plot(x, Y)
points(Y_hat$x, Y_hat$y_hat, type="l")

## KDE
p_hat<-p(Y, h=1/5)
plot(density(Y), col="orange")
points(p_hat$x, p_hat$p_hat, type="l", col="blue")

# Poisson Example
## Non-parametric regression
n<-100
x<-seq(-2,2,length.out=n)
X<-matrix(c(rep(1,n), x),ncol=2)
b<-c(.922,-.501)
lam<-exp(X%*%b)
Y<-c()
for(i in 1:n){
  Y<-c(Y, rpois(1, lam[i]))
}

Y_hat<-r(x, Y)
plot(x, Y)
points(Y_hat$x, Y_hat$y_hat, type="l")

# KDE
p_hat<-p(Y, h=1/2)
plot(density(Y), col="orange")
points(p_hat$x, p_hat$p_hat, type="l", col="blue")


# Alcohol Example
## Read data
mvalc <- read.csv("/home/pat-mellady/Desktop/Tidy gark/alc_data.csv", header=T)

## Tidy Data
rownames(mvalc) <- mvalc$Wedding
mvalc<-t(mvalc[,-1])

## Define Variables
X<-mvalc[,1]
X2<-X^2

Y<-mvalc[,62]
Y_w<-mvalc[,63]

high<-mvalc[c(2,6,7,9,12,13),c(1,62)]
low<-mvalc[-c(2,6,7,9,12,13), c(1,62)]

Y_h<-high[,2]
X_h<-high[,1]
X_h2<-X_h^2

Y_l<-low[,2]
X_l<-low[,1]
X_l2<-X_l^2

## Non-parametric regression
x_h<-seq(min(X_h), max(X_h), length.out=10*length(X_h))
Y_h_hat<-r(X_h, Y_h, h=10)
plot(X_h, Y_h, col="blue")
points(Y_h_hat$x, Y_h_hat$y_hat, type="l", col="blue")

x_l<-seq(min(X_l), max(X_l), length.out=10*length(X_l))
Y_l_hat<-r(X_l, Y_l, h=1)
points(X_l, Y_l, col="orange")
points(Y_l_hat$x, Y_l_hat$y_hat, type="l", col="orange")
