# Set-up -----------------------------------------------------------------------
library(missForest)
library(MASS)
library(mgcv)
library(pg)

# Necessary transformation functions -------------------------------------------
## Function to transform from psi back to pi
stick_breaking_pi<-function(psi){
  Q<-length(psi)
  pi<-numeric(Q+1)
  remaining<-1
  for(k in 1:Q){
    tilde_pi_k<-1/(1+exp(-psi[k]))   # sigma(psi_k)
    pi[k]<-tilde_pi_k*remaining
    remaining<-remaining-pi[k]
  }
  pi[Q+1]<-remaining
  pi
}

## Function to calculate the probability of interest
estimate_pi<-function(data, samples, burnin){
  n<-nrow(data$Y)
  nsamps<-length(samples$B)
  burn_idx<-(burnin+1):nsamps
  
  X<-data$X_fixed
  Z<-data$X_random
  
  Q<-ncol(data$Y)
  p_x<-ncol(X)
  q<-ncol(Z)
  
  pi_hat<-vector("list", length(burn_idx))
  p_hat<-matrix(0, nrow=n, ncol=length(burn_idx))
  idx<-0
  for(r in burn_idx){
    idx<-idx+1
    
    B_mat<-matrix(samples$B[[r]], nrow=p_x, ncol=Q)
    b_mat<-matrix(samples$b[[r]], nrow=q, ncol=Q)
    psi_hat_r<-X%*%B_mat+Z%*%b_mat
    
    pi_hat_r<-t(apply(psi_hat_r, 1, stick_breaking_pi))
    p_hat_r<-pi_hat_r[,1]/(pi_hat_r[,1]+pi_hat_r[,4])
    
    pi_hat[[idx]]<-pi_hat_r
    p_hat[,idx]<-p_hat_r
    
  }
  
  p_hat_mean<-rowMeans(p_hat)
  return(list(point_est=p_hat_mean, p_hat=p_hat, pi_hat=pi_hat))
  
}

# Helper functions for data/sample containers ----------------------------------
## Function to help format the data nicely for use in the sampler
build_matrices_1<-function(df, K, fixed_vars, random_vars){
  fixed_formula_text<-paste0("~ ", paste0(fixed_vars, collapse=" + "))
  fixed_formula<-as.formula(fixed_formula_text)
  X_fixed<-model.matrix(fixed_formula, df)
  
  random_formula_text<-paste0("~ ", paste0(random_vars, collapse=" + "))
  random_formula<-as.formula(random_formula_text)
  X_random<-model.matrix(random_formula, df)
  
  n<-nrow(df)
  Y<-matrix(0, nrow=n, ncol=K-1)
  for(i in 1:n){
    if(df[i,"is_strike"]==1 & df[i,"is_swing"]==0){
      Y[i,1]<-1
    } else if(df[i,"is_strike"]==1 & df[i,"is_swing"]==1){
      Y[i,2]<-1
    } else if(df[i,"is_strike"]==0 & df[i,"is_swing"]==1){
      Y[i,3]<-1
    }
  }
  
  N<-matrix(0, nrow=n, ncol=K-1)
  N[,1]<-1
  if(K-1>=2){
    for(k in 2:(K-1)){
      N[,k]<-N[,1]-rowSums(Y[,1:(k-1),drop=FALSE])
    }
  }
  
  return(list(Y=Y, N=N, X_fixed=X_fixed, X_random=X_random))
}

## Function to make a list-of-lists style sample container
initialize_samples<-function(nsamps, n, Q, p_x, q){
  B<-vector("list", nsamps)
  for(r in 1:nsamps){
    B[[r]]<-vector("numeric", Q*p_x)
    if(r==1){
      B[[r]]<-rep(0, Q*p_x)
    }
  }
  
  b<-vector("list", nsamps)
  for(r in 1:nsamps){
    b[[r]]<-vector("numeric", Q*q)
    if(r==1){
      b[[r]]<-rep(0, Q*q)
    }
  }
  
  omega<-vector("list", n)
  for(i in 1:n){
    omega[[i]]<-rep(1, Q)
  }
  
  return(list(B=B, b=b, omega=omega))
}

# Main sampler -----------------------------------------------------------------
## Log-likelihood function to provide diagnostic prints
log_lik<-function(data, samples, r){
  Y<-data$Y
  X<-data$X_fixed
  Z<-data$X_random
  
  Q<-ncol(Y)
  p_x<-ncol(X)
  q<-ncol(Z)
  
  B_r<-matrix(samples$B[[r]], nrow=p_x, ncol=Q)
  b_r<-matrix(samples$b[[r]], nrow=q, ncol=Q)
  
  
  psi_r<-X%*%B_r+Z%*%b_r
  pi_r<-t(apply(psi_r, 1, stick_breaking_pi))
  
  Y_aug<-cbind(Y, 1-rowSums(Y))
  
  log_pi_r<-log(pi_r+1e-16)
  
  temp<-Y_aug*log_pi_r
  ll<-sum(temp)
  
  return(ll)
}

## Brier score function to calculate and track the final metric
brier_score<-function(df, data, samples, r){
  n<-nrow(data$Y)
  
  X<-data$X_fixed
  Z<-data$X_random
  
  Q<-ncol(data$Y)
  p_x<-ncol(X)
  q<-ncol(Z)
  
  B_mat<-matrix(samples$B[[r]], nrow=p_x, ncol=Q)
  b_mat<-matrix(samples$b[[r]], nrow=q, ncol=Q)
  psi_hat_r<-X%*%B_mat+Z%*%b_mat
  
  pi_hat_r<-t(apply(psi_hat_r, 1, stick_breaking_pi))
  p_hat<-pi_hat_r[,1]/(pi_hat_r[,1]+pi_hat_r[,4])
  
  brier<-mean((p_hat-df[["is_event"]])^2, na.rm=TRUE)
  
  return(brier)
}

## Original sampler code (slow)
MN_GIBB<-function(nsamps, data, B_0, Sigma_B, b_0, Sigma_b, trace=TRUE){
  
  ## Extract necessary constants
  mu_r<-data$Y-data$N/2
  Q<-ncol(mu_r)
  p_x<-ncol(data$X_fixed)
  q<-ncol(data$X_random)
  n<-nrow(data$Y)
  
  ## Pull design matrices
  X<-data$X_fixed
  Z<-data$X_random
  
  ## Initialize samples
  samples<-initialize_samples(nsamps, n, Q, p_x, q)
  
  ## Initialize psi
  psi_r<-matrix(0, nrow=n, ncol=Q)
  ll<-c(log_lik(data, samples, 1))
  bs<-c(brier_score(df, data, samples, 1))
  msg_interval<-max(1, floor(nsamps/10))
  # Begin sampling
  if(trace) cat("Running sampler...\n")
  for(r in 2:nsamps){
    
    # Update mu
    for(i in 1:n){
      mu_r[i,]<-(data$Y[i,]-data$N[i,]/2)/samples$omega[[i]]
    }
    
    # Sample the fixed effect coefficients
    temp<-rep(0, Q*p_x)
    Sigma_tilde_inv<-matrix(0, ncol=Q*p_x, nrow=Q*p_x)
    for(i in 1:n){
      temp_kron_X<-t(kronecker(diag(Q), X[i,]))
      temp_kron_Z<-t(kronecker(diag(Q), Z[i,]))
      Omega_rk<-diag(drop(samples$omega[[i]]), ncol=Q, nrow=Q)
      Sigma_tilde_inv<-Sigma_tilde_inv+t(temp_kron_X)%*%Omega_rk%*%temp_kron_X
      temp<-temp+t(temp_kron_X)%*%Omega_rk%*%(mu_r[i,]-temp_kron_Z%*%samples$b[[r-1]])
    }
    Sigma_tilde<-solve(Sigma_tilde_inv+solve(Sigma_B))
    B_tilde<-Sigma_tilde%*%(temp+solve(Sigma_B)%*%B_0)
    samples$B[[r]]<-mvrnorm(1, B_tilde, Sigma_tilde)
    
    # Sample the random effect coefficients
    temp<-rep(0, Q*q)
    Sigma_tilde_inv<-matrix(0, ncol=Q*q, nrow=Q*q)
    for(i in 1:n){
      temp_kron_X<-t(kronecker(diag(Q), X[i,]))
      temp_kron_Z<-t(kronecker(diag(Q), Z[i,]))
      Omega_rk<-diag(drop(samples$omega[[i]]), ncol=Q, nrow=Q)
      Sigma_tilde_inv<-Sigma_tilde_inv+t(temp_kron_Z)%*%Omega_rk%*%temp_kron_Z
      temp<-temp+t(temp_kron_Z)%*%Omega_rk%*%(mu_r[i,]-temp_kron_X%*%samples$B[[r]])
    }
    Sigma_tilde<-solve(Sigma_tilde_inv+solve(Sigma_b))
    b_tilde<-Sigma_tilde%*%(temp+solve(Sigma_b)%*%b_0)
    samples$b[[r]]<-mvrnorm(1, b_tilde, Sigma_tilde)
    
    # Update psi
    for(i in 1:n){
      temp_kron_X<-kronecker(diag(Q), X[i,])
      temp_kron_Z<-kronecker(diag(Q), Z[i,])
      psi_r[i,]<-temp_kron_X%*%samples$B[[r]]+temp_kron_Z%*%samples$b[[r]]
    }
    
    # Sample the augmenting PG variables
    for(i in 1:n){
      for(k in 1:Q){
        samples$omega[[i]][k]<-rpg_scalar(data$N[i,k], psi_r[i,k])
      }
    }
    
    ## Track log likelihood and Brier score
    ll<-c(ll, log_lik(data, samples, r))
    bs<-c(bs, brier_score(df, data, samples, r))
    
    ## Print progress
    if(trace){
      if(r%%msg_interval==0 | r==nsamps){
        pct_complete<-round((r/nsamps)*100)
        ll_r<-ll[r]
        bs_r<-bs[r]
        bar_width<-20
        filled<- floor((r/nsamps)*bar_width)
        bar<-paste0("|", strrep("=", filled), strrep("-", bar_width-filled), "|")
        
        cat(sprintf("\r%3d%% %s Iter: %d/%d | LL: %s | BS: %s",
                    pct_complete,
                    bar,
                    r,
                    nsamps,
                    format(round(ll_r, 1), big.mark=","),
                    format(round(bs_r, 3), big.mark=",")))
        
        utils::flush.console()
        
        if(r==nsamps) cat("\nDone.\n")
      }
    }
  }
  
  return(samples)
}

## Speed boosted sampler code
MN_GIBB_2<-function(nsamps, df, data, B_0, Sigma_B, b_0, Sigma_b, trace=TRUE){
  ## Extract necessary constants
  mu_r<-data$Y-data$N/2
  Q<-ncol(mu_r)
  p_x<-ncol(data$X_fixed)
  q<-ncol(data$X_random)
  n<-nrow(data$Y)
  
  ## Pull design matrices
  X<-data$X_fixed
  Z<-data$X_random
  
  ## Initialize samples
  samples<-initialize_samples(nsamps, n, Q, p_x, q)
  
  ## Initialize psi
  psi_r<-matrix(0, nrow=n, ncol=Q)
  ll<-c(log_lik(data, samples, 1))
  bs<-c(brier_score(df, data, samples, 1))
  msg_interval<-max(1, floor(nsamps/10))
  # Begin sampling
  if(trace) cat("Running sampler...\n")
  for(r in 2:nsamps){
    
    ## Aggregate all current omega values into one matrix
    Om<-do.call(rbind, samples$omega)
    ## Calculate kappa (tranformed mu)
    kappa<-data$Y-data$N/2
    
    ## Pull previous iteration values into storage for easier computation
    B_prev<-matrix(samples$B[[r-1]], nrow=p_x, ncol=Q)
    b_prev<-matrix(samples$b[[r-1]], nrow=q, ncol=Q)
    B_new<-matrix(0, p_x, Q)
    b_new<-matrix(0, q, Q)
    
    for(k in 1:Q){
      ## Get the kth omega vector
      om_k<-Om[,k]
      
      ## Update the fixed effects
      XtOX<-crossprod(X, X*om_k)
      XtOZb<-crossprod(X, om_k*(Z%*%b_prev[,k]))
      Xtkappa<-crossprod(X, kappa[,k])
      
      Sig_B_inv_k<-solve(Sigma_B[((k-1)*p_x+1):(k*p_x), ((k-1)*p_x+1):(k*p_x)])
      Sigma_tilde<-solve(XtOX+Sig_B_inv_k)
      B_tilde<-Sigma_tilde%*%(Xtkappa-XtOZb+Sig_B_inv_k%*%B_0[((k-1)*p_x+1):(k*p_x)])
      B_new[,k]<-mvrnorm(1, B_tilde, Sigma_tilde)
      
      ## Update the random effects
      ZtOZ<-crossprod(Z, Z*om_k)
      ZtOXB<-crossprod(Z, om_k*(X%*%B_new[,k]))
      Ztkappa<-crossprod(Z, kappa[,k])
      
      Sig_b_inv_k<-solve(Sigma_b[((k-1)*q+1):(k*q), ((k-1)*q+1):(k*q)])
      Sigma_tilde_b<-solve(ZtOZ+Sig_b_inv_k)
      b_tilde<-Sigma_tilde_b%*%(Ztkappa-ZtOXB+Sig_b_inv_k%*%b_0[((k-1)*q+1):(k*q)])
      b_new[,k]<-mvrnorm(1, b_tilde, Sigma_tilde_b)
    }
    
    ## Store the updated parameters in the list
    samples$B[[r]]<- as.vector(B_new)
    samples$b[[r]]<- as.vector(b_new)
    
    ## Calculate psi and sample new omega
    psi<-X%*%B_new+Z%*%b_new
    omega_new<-vector("list", n)
    for(i in 1:n){
      omega_new[[i]]<-vector("numeric", Q)
      for(k in 1:Q){
        omega_new[[i]][k]<-rpg_scalar(data$N[i,k], psi[i,k])
      }
    }
    samples$omega<-omega_new
    
    ## Track log likelihood and Brier score
    ll<-c(ll, log_lik(data, samples, r))
    bs<-c(bs, brier_score(df, data, samples, r))
    
    ## Print progress
    if(trace){
      if(r%%msg_interval==0 | r==nsamps){
        pct_complete<-round((r/nsamps)*100)
        ll_r<-ll[r]
        bs_r<-bs[r]
        bar_width<-20
        filled<- floor((r/nsamps)*bar_width)
        bar<-paste0("|", strrep("=", filled), strrep("-", bar_width-filled), "|")
        
        cat(sprintf("\r%3d%% %s Iter: %d/%d | LL: %s | BS: %s",
                    pct_complete,
                    bar,
                    r,
                    nsamps,
                    format(round(ll_r, 1), big.mark=","),
                    format(round(bs_r, 3), big.mark=",")))
        
        utils::flush.console()
        
        if(r==nsamps) cat("\nDone.\n")
      }
    }
  }
  
  return(list(log_lik=ll, bs=bs, samples=samples))
}



