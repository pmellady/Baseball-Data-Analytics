################################################################################
########## Neural Network PDE Solver for Burgers' Equation #####################
################################################################################

# Setting Parameters
n<-32                                  # n^2 is the number of meshpoints in the region
d<-2                                    # the number of dimensions of the input
q<-10                                    # the number of neurons in the hidden layer

# Initialize weights on the nodes
U<-matrix(rep(1, q), ncol=1)            
V<-matrix(rep(1, q), ncol=1)
W<-matrix(rep(1,q*d), ncol=d)

# Components of the network ####################################################
# Sigmoid function for the activation of the neurons
s<-function(x){
  res<-1/(1+exp(-x))
  return(res)
}

# Derivative of the sigmoid function
# (used in derivatives of the estimated solution)
s_prime<-function(x){
  res<-exp(-x)/(1+exp(-x))^2
  return(res)
}

# Second derivative of the sigmoid function
# (used in derivatives of the estimated solution)
s_dprime<-function(X){
  n<-length(X)
  res<-c()
  for(i in 1:n){
    temp1<--exp(-X[i])/(1+exp(-X[i]))^2
    temp2<-2*exp(-2*X[i])/(1+exp(-X[i]))^3
    temp<-temp1+temp2
    res<-c(res,temp)
  }
  return(matrix(res, ncol=1))
}

N<-function(X, U, V, W){
  N_x<-t(V)%*%s(U+W%*%X)
  return(N_x)
  
}

# Partials of the non-unique portion of the estimated solution
N_t<-function(X, U, V, W){
  N_t_xt<-sum(V*W[,2]*s_prime(U+W%*%X))
  return(N_t_xt)
}

N_x<-function(X, U, V, W){
  N_x_xt<-sum(V*W[,1]*s_prime(U+W%*%X))
  return(N_x_xt)
}

E<-function(U, V, W){
  grid<-seq(0,1,length.out=n)
  temp<-c()
  for(x in grid){
    for(t in grid){
      X<-matrix(c(x,t))
      temp<-c(temp, (N(X, U, V, W)+x*N_x(X, U, V, W)+(x*N(X, U, V, W)+t^3)*(x*N_t(X, U, V, W)+3*t^2))^2)
    }
  }
  return(sum(temp))
}

g_mu<-function(U, V, W, du, dv, dw, mu){
  gmu_u<-(1/mu)*(E(U+mu*du, V, W)-E(U, V, W))[[1]]*du
  gmu_v<-(1/mu)*(E(U, V+mu*dv, W)-E(U, V, W))[[1]]*dv
  gmu_w<-(1/mu)*(E(U, V, W+mu*dw)-E(U, V, W))[[1]]*dw
  return(list(U_k=gmu_u, V_k=gmu_v, W_k=gmu_w))
}

# Back tracking line search to find optimal step size
back_track<-function(U, V, W, g_k, rho=0.3, c=0.005, h=4.05e-6){
  t<-0.005
  while(E(U-h*g_k[[1]], V-h*g_k[[2]], W-h*g_k[[3]])-E(U, V, W)>=c*t){
    h<-h*rho
  }
  return(h)
}

sol<-function(U0, V0, W0, K=5000, h0=4.05e-6, mu=4.05e-8, epsilon=0.5){
  k<-0
  error<-E(U0, V0, W0)
  cat("0 garkages with a finna garkage of ", error,"\n")
  while(error>epsilon){
    du<-matrix(rnorm(q, 0, 1), ncol=1)
    dv<-matrix(rnorm(q, 0, 1), ncol=1)
    dw<-matrix(rnorm(d*q, 0, 1), ncol=d)
    g_k<-g_mu(U0, V0, W0, du, dv, dw, mu)
    h<-back_track(U0, V0, W0, g_k, h=h0)
    U_k<-U0-h*g_k[[1]]
    V_k<-V0-h*g_k[[2]]
    W_k<-W0-h*g_k[[3]]
    U0<-U_k
    V0<-V_k
    W0<-W_k
    error<-E(U_k, V_k, W_k)
    k<-k+1
    if(k%%25==0){
      cat(k, " garkages with a finna garkage of ", error, " and a stepping gark of ", h ,"\n")
    }
    if(k>K){
      break
    }
  }
  return(list(U=U_k, V=V_k, W=W_k))
}

# Warm start to speed up the convergence
sol_warm_start<-function(U0, V0, W0, K=1000, epsilon=0.8, mu=4.05e-8, h0=4.05e-8, warm=T, warm_iter=300){
  i<-0
  if(warm){
    cat("Warm finnage","\n")
    u<-U0
    v<-V0
    w<-W0
    for (i in 1:warm_iter) {
      p<-.4
      variance <- .2
      stepu<-p*matrix(rnorm(q, 0, variance), ncol=1)
      stepv<-p*matrix(rnorm(q, 0, variance), ncol=1)
      stepw<-p*matrix(rnorm(d*q, 0, variance), ncol=d)
      
      unew <- u + stepu
      wnew <- w + stepw
      vnew <- v + stepv
      counter <- 0
      
      unew1 <- unew
      wnew1 <- wnew
      vnew1 <- vnew
      
      while(E(unew1, vnew1, wnew1) < E(u, v, w)) {
        unew1 <- unew1 + stepu
        wnew1 <- wnew1 + stepw
        vnew1 <- vnew1 + stepv
        counter <- counter + 1
      }
      unew <- unew1 - stepu
      wnew <- wnew1 - stepw
      vnew <- vnew1 - stepv
      if (E(unew, vnew, wnew) < E(u, v, w)) {
        u <- unew
        v <- vnew
        w <- wnew
      }
      cat("Level ", counter, " finna garkage!!", E(u, v, w), " \n")
      print(i)
    }
    U0<-u
    V0<-v
    W0<-w
  }
  k<-0
  error<-E(U0, V0, W0)
  cat("0 garkages with a finna garkage of ", error,"\n")
  while(error>epsilon){
    du<-matrix(rnorm(q, 0, 1), ncol=1)
    dv<-matrix(rnorm(q, 0, 1), ncol=1)
    dw<-matrix(rnorm(d*q, 0, 1), ncol=d)
    g_k<-g_mu(U0, V0, W0, du, dv, dw, mu)
    h<-back_track(U0, V0, W0, g_k, h=h0)
    U_k<-U0-h*g_k[[1]]
    V_k<-V0-h*g_k[[2]]
    W_k<-W0-h*g_k[[3]]
    U0<-U_k
    V0<-V_k
    W0<-W_k
    error<-E(U_k, V_k, W_k)
    k<-k+1
    if(k%%25==0){
      cat(k, " Nesterov garkages with a finna garkage of ", error, " and a stepping gark of ", h ,"\n")
    }
    if(k>K){
      break
    }
  }
  return(list(U=U0, V=V0, W=W0))
}


# Let's solve some PDE's #######################################################
params<-sol_warm_start(U, V, W, epsilon = 0.01, warm=T, h0=.05)
U<-params[[1]]
V<-params[[2]]
W<-params[[3]]
params<-sol_warm_start(U, V, W, K=5000, epsilon = 0.1, warm=F, h0=1, mu=4.05e-3)
U<-params[[1]]
V<-params[[2]]
W<-params[[3]]


# Plotting #####################################################################
U0<-matrix(rep(1, q), ncol=1)            
V0<-matrix(rep(1, q), ncol=1)
W0<-matrix(rep(1,q*d), ncol=d)


grid<-seq(0,1,length.out=n)
#region<-matrix(rep(0, n^2), ncol=n)
#Y<-matrix(rep(0, n^2), ncol=n)
region<-c()
Y<-c()
og<-c()
for(i in 1:n){
  temp1<-c()
  temp2<-c()
  temp3<-c()
  for(j in 1:n){
    x<-grid[i]
    t<-grid[j]
    #region[i,j]<-x*N(c(x, t), U, V, W)+t^3
    #Y[i,j]<-exp(-3*x)*t^3
    temp1<-c(temp1,x*N(c(x, t), U, V, W)+t^3)
    temp2<-c(temp2,exp(-3*x)*t^3)
    temp3<-c(temp3, x*N(c(x, t), U0, V0, W0)+t^3)
  }
  region<-c(region, temp1)
  Y<-c(Y, temp2)
  og<-c(og, temp3)
  
}

dim(region)<-c(32,32)
dim(Y)<-c(32,32)
dim(og)<-c(32,32)

persp(grid, grid, region, xlab="x", ylab="t", zlab="y_hat")
persp(grid, grid, Y, xlab="x", ylab="t", zlab="y")



U<-c(24.575445, -93.709455, 51.992195, 13.740927, -34.342413, -23.587034, 81.668723, -23.444581, -6.297116, 48.975537)
V<-c(1.067274, -94.545727, -12.115458, 7.659357, 5.128956, 45.461944, 12.869272, -23.557857, 8.306440, -9.617682)
W<-matrix(c(8.456339, -33.886397, -27.424763,  15.086468, 85.641890, -45.017892, 9.931575,  25.425734, -38.464277,  -8.331545,
            -2.052845, -56.874900, 14.323503,  31.641610, -30.667624, -11.181155, -19.980516, -29.612373, -20.221208,  19.574767), ncol=2, byrow=TRUE)


library(plotly)
fig <- plot_ly(showscale = FALSE)
fig <- fig %>% add_surface(z = ~Y)
fig <- fig %>% add_surface(z = ~region, opacity = 0.98)


fig
