################################################################################
########## Neural Network ODE Solver for the equation y''+2y'+y=2e^(-x) ########
################################################################################

# Setting Parameters
n<-1000                                     # the number of mesh points in the interval
q<-10                                       # number of neurons (one layer)
X<-matrix(seq(0,1,length.out=n), ncol=1)    # making the grid of mesh points

# Initialize weights on the nodes
U<-matrix(rep(1, q), ncol=1)
V<-matrix(rep(1, q), ncol=1)
W<-matrix(rep(1,q), ncol=1)

# Components of the network ####################################################
# Sigmoid function for the activation of the neurons
s<-function(X){
  n<-length(X)
  res<-c()
  for(i in 1:n){
    temp<-1/(1+exp(-X[i]))
    res<-c(res, temp)
  }
  return(res)
}

# Derivative of the sigmoid function
# (used in derivatives of the estimated solution)
s_prime<-function(X){
  n<-length(X)
  res<-c()
  for(i in 1:n){
    temp<-exp(-X[i])/(1+exp(-X[i]))^2
    res<-c(res,temp)
  }
  return(matrix(res, ncol=1))
}

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

# Non-unique portion of the estimated solution #################################
# Defining the non-unique portion of the estimated solution
# (does not include the boundary conditions)
N<-function(x, U, V, W){
  S<-s(U+W*x)
  N_x<-t(V)%*%S
  return(N_x)
  
}

# Defining the derivative of the non-unique portion of the estimated solution
# (does not include the boundary conditions)
N_prime<-function(x, U, V, W){
  N_prime_x<-t(W)%*%(V*s_prime(U+W*x))
  return(N_prime_x)
}

N_dprime<-function(x, U, V, W){
  N_dprime_x<-sum(V*W^2*s_dprime(U+W*x))
  return(N_dprime_x)
}

# Changing the problem to an optimization problem ##############################
# Defining the error function that we want to minimize (includes the boundary conditions)
E<-function(X, U, V, W){
  n<-length(X)
  temp<-c()
  for(i in 1:n){
    x<-X[i]
    temp<-c(temp, ((x^2+4*x+2)*N(x, U, V, W)+(2*x^2+4*x)*N_prime(x, U, V, W)+x^2*N_dprime(x, U, V, W)-2*x*exp(x)+22*x+53)^2)
  }
  return(sum(temp))
}


# Components of the gradient-free method #######################################
# Calculating the numerical directional derivative for the optimization method
g_mu<-function(x, U_k, V_k, W_k, du, dv, dw, mu){
  gmu_u<-(1/mu)*(E(x, U_k+mu*du, V_k, W_k)-E(x, U_k, V_k, W_k))[[1]]*du
  gmu_v<-(1/mu)*(E(x, U_k, V_k+mu*dv, W_k)-E(x, U_k, V_k, W_k))[[1]]*dv
  gmu_w<-(1/mu)*(E(x, U_k, V_k, W_k+mu*dw)-E(x, U_k, V_k, W_k))[[1]]*dw
  return(list(U_k=gmu_u, V_k=gmu_v, W_k=gmu_w))
}

# Back tracking line search to find optimal step size
back_track<-function(X, U_k, V_k, W_k, g_k, rho=0.3, c=0.005, h=4.05e-6){
  t<-0.005
  while(E(X, U_k-h*g_k[[1]], V_k-h*g_k[[2]], W_k-h*g_k[[3]])-E(X, U_k, V_k, W_k)>=c*t){
    h<-h*rho
  }
  return(h)
}

U0<-U
V0<-V
W0<-W

sol<-function(X, U0, V0, W0, h0=4.05e-6, mu=4.05e-8, epsilon=0.5){
  K<-5000
  k<-0
  error<-E(X, U0, V0, W0)
  cat("0 garkages with a finna garkage of ", error,"\n")
  while(error>epsilon){
    du<-matrix(rnorm(q, 0, 1), ncol=1)
    dv<-matrix(rnorm(q, 0, 1), ncol=1)
    dw<-matrix(rnorm(q, 0, 1), ncol=1)
    g_k<-g_mu(X, U0, V0, W0, du, dv, dw, mu)
    h<-back_track(X, U0, V0, W0, g_k, h=h0)
    U_k<-U0-h*g_k[[1]]
    V_k<-V0-h*g_k[[2]]
    W_k<-W0-h*g_k[[3]]
    U0<-U_k
    V0<-V_k
    W0<-W_k
    error<-E(X, U_k, V_k, W_k)
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
sol_warm_start<-function(X, U0, V0, W0, K=1000, epsilon=0.8, mu=4.05e-8, h=4.05e-8, warm=T, warm_iter=150){
  n<-length(X)
  i<-0
  if(warm){
    u<-U0
    v<-V0
    w<-W0
    for (i in 1:warm_iter) {
      p<-0.1
      variance <- .2
      stepu<-p*matrix(rnorm(q, 0, variance), ncol=1)
      stepv<-p*matrix(rnorm(q, 0, variance), ncol=1)
      stepw<-p*matrix(rnorm(q, 0, variance), ncol=1)
      
      unew <- u + stepu
      wnew <- w + stepw
      vnew <- v + stepv
      counter <- 0
      
      unew1 <- unew
      wnew1 <- wnew
      vnew1 <- vnew
      
      while(E(X, unew1, vnew1, wnew1) < E(X, u, v, w)) {
        unew1 <- unew1 + stepu
        wnew1 <- wnew1 + stepw
        vnew1 <- vnew1 + stepv
        counter <- counter + 1
      }
      unew <- unew1 - stepu
      wnew <- wnew1 - stepw
      vnew <- vnew1 - stepv
      if (E(X, unew, vnew, wnew) < E(X, u, v, w)) {
        u <- unew
        v <- vnew
        w <- wnew
      }
      cat("Level ", counter, " finna garkage!!", E(X, u, v, w), " \n")
      print(i)
    }
    U0<-u
    V0<-v
    W0<-w
  }
  error<-E(X, U0, V0, W0)
  k<-0
  while(error>epsilon){
    du<-matrix(rnorm(q, 0, 1), ncol=1)
    dv<-matrix(rnorm(q, 0, 1), ncol=1)
    dw<-matrix(rnorm(q, 0, 1), ncol=1)
    g_k<-g_mu(X, U0, V0, W0, du, dv, dw, mu)
    h<-back_track(X, U0, V0, W0, g_k)
    U_k<-U0-h*g_k[[1]]
    V_k<-V0-h*g_k[[2]]
    W_k<-W0-h*g_k[[3]]
    U0<-U_k
    V0<-V_k
    W0<-W_k
    error<-E(X, U_k, V_k, W_k)
    k<-k+1
    if(k%%25==0){
      cat(k, " Nesterov garkages with a finna garkage of ", error, " and a stepping gark of ", h ,"\n")
    }
    if(k>5000){
      break
    }
  }
  return(list(U=U_k, V=V_k, W=W_k))
}


# Let's solve some ODE's #######################################################
params<-sol_warm_start(X, U, V, W, epsilon = 1000)
U<-params[[1]]
V<-params[[2]]
W<-params[[3]]
params<-sol(X, U, V, W, epsilon = 1000)
U<-params[[1]]
V<-params[[2]]
W<-params[[3]]
params<-sol(X, U, V, W, epsilon = 100, h0=0.005)
U<-params[[1]]
V<-params[[2]]
W<-params[[3]]
params<-sol(X, U, V, W, epsilon = 25, h0=0.05)

y_hat<-c()
for(i in 1:n){
  temp<-X[i]^2*N(X[i], U, V, W)+X[i]*22+9
  y_hat<-c(y_hat, temp)
  
}
curve(9*exp(-x)+31*x*exp(-x)+x^3*exp(-x)/8, col="blue")
points(X, y_hat, type="l", col="orange")




U<-c(27.900562, 4.668395, 86.398555, -120.130045, 16.396949, -67.857932, -145.844922, 194.924957, 18.869293, -1.771952)
V<-c(36.32949, 137.90982, -183.43228, -139.29383, 50.14412, 68.43253, 132.70959, -211.60481, 131.93258, 111.84776)
W<-c(146.6507617, 66.2037910, -63.7876460, 10.3580539, 32.0746499, -101.7057346, -110.9633051,-34.7742198, 211.4974973, 0.4143304)
