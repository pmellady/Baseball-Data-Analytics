# Set-up -----------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(data.table)

# Helper functions to perform optimization -------------------------------------
## Define function to transform linear predictor to estimate pi
calc_pi<-function(n, m, Y, S, X, rho, B, K, autoregress=TRUE, exogenous=TRUE){
  J<-nrow(Y[[1]])
  B_r<-matrix(B, ncol=K)
  pi<-vector("list", n)
  for(i in 1:n){
    pi[[i]]<-matrix(0, nrow=J, ncol=K)
    for(j in 1:J){
      X_ij<-X[[i]][j,]
      if(exogenous){
        eta<-X_ij%*%B_r
        if(autoregress & j>1){
          eta<-eta+as.numeric(rho*S[[i]][j-1,])
        }
      } else{
        if(autoregress){
          if(j==1){
            next
          } else{
            eta<-as.numeric(rho%*%Y[[i]][j-1,]/m[i,j-1])
          }
        } else{
          eta<-matrix(0, ncol=K, nrow=J)
        }
      }
      pi[[i]][j,]<-exp(eta)/(1+sum(exp(eta)))
    }
  }
  return(pi)
}

## Function to calculate gradient of log-likelihood w.r.t B
nabla_B<-function(n, X, Y, m, B, pi_hat){
  k<-ncol(Y[[1]])
  J<-nrow(Y[[1]])
  
  grad<-numeric(length(B))
  for(i in 1:n){
    for(j in 1:J){
      resid<-as.matrix(Y[[i]][j,]-m[i,j]*pi_hat[[i]][j,])
      grad<-grad+kronecker(diag(k), X[[i]][j,])%*%resid
    }
  }
  return(grad)
}

## Function to calculate the Hessian of log-likelihood w.r.t B
nabla_2_B<-function(n, X, Y, m, B, pi_hat){
  J<-nrow(Y[[1]])
  k<-ncol(Y[[1]])
  
  hess<-matrix(0, ncol=length(B), nrow=length(B))
  for(i in 1:n){
    for(j in 1:J){
      resid<-Y[[i]][j,]-m[i,j]*pi_hat[[i]][j,]
      grad<-kronecker(diag(k), X[[i]][j,])%*%resid
      hess<-hess-grad%*%t(grad)
    }
  }
  return(hess)
}

## Function to calculate gradient of log-likelihood w.r.t rho
nabla_rho<-function(n, S, Y, m, pi_hat){
  J<-nrow(Y[[1]])
  K<-ncol(Y[[1]])
  grad<-rep(0, K)
  for(i in 1:n){
    for(j in 2:J){
      resid<-as.matrix(Y[[i]][j,]-m[i,j]*pi_hat[[i]][j,])
      scale<-S[[i]][j-1,]
      grad<-grad+scale*resid
    }
  }
  return(grad)
}

## Function to calculate the Hessian of log-likelihood w.r.t rho
nabla_2_rho<-function(n, S, Y, m, pi_hat){
  J<-nrow(Y[[1]])
  K<-ncol(Y[[1]])
  hess<-matrix(0, ncol=K, nrow=K)
  for(i in 1:n){
    for(j in 2:J){
      resid<-as.matrix(Y[[i]][j,]-m[i,j]*pi_hat[[i]][j,])
      scale<-S[[i]][j-1,]
      grad_j<-scale*resid
      hess<-hess-grad_j%*%t(grad_j)
    }
  }
  return(hess)
}

## Function to calculate gradient of log-likelihood w.r.t rho
nabla_rho_matrix<-function(n, S, Y, m, pi_hat){
  J<-nrow(Y[[1]])
  K<-ncol(Y[[1]])
  grad<-rep(0, K^2)
  for(i in 1:n){
    for(j in 2:J){
      resid<-as.matrix(Y[[i]][j,]-m[i,j]*pi_hat[[i]][j,])
      scale<-kronecker((Y[[i]][j-1,]/m[i,j-1]), diag(K))
      grad<-grad+scale%*%resid
    }
  }
  return(grad)
}

## Function to calculate the Hessian of log-likelihood w.r.t rho
nabla_2_rho_matrix<-function(n, S, Y, m, pi_hat){
  J<-nrow(Y[[1]])
  K<-ncol(Y[[1]])
  hess<-matrix(0, ncol=K^2, nrow=K^2)
  for(i in 1:n){
    for(j in 2:J){
      resid<-as.matrix(Y[[i]][j,]-m[i,j]*pi_hat[[i]][j,])
      scale<-kronecker((Y[[i]][j-1,]/m[i,j-1]), diag(K))
      grad_j<-scale%*%resid
      hess<-hess-grad_j%*%t(grad_j)
    }
  }
  return(hess)
}

# Helper function to format the data nicely ------------------------------------
## Function to build the data used in fitting
build_matrices<-function(df, subject_id, y_cols, s_cols, fixed_vars){
  ids<-sort(unique(df[[subject_id]]))
  n<-length(ids)
  
  fixed_formula_text<-paste0("~ ", paste0(fixed_vars, collapse=" + "))
  fixed_formula<-as.formula(fixed_formula_text)
  
  X<-vector("list", n)
  Y<-vector("list", n)
  S<-vector("list", n)
  m<-NULL
  
  for(i in seq_along(ids)){
    id<-ids[i]
    subj_df<-df[df[[subject_id]]==id, ]
    
    m<-rbind(m, subj_df[["m"]])
    X[[i]]<-model.matrix(fixed_formula, data=subj_df)
    S[[i]]<-as.matrix(subj_df[,s_cols])
    Y[[i]]<-as.matrix(subj_df[,y_cols])
  }
  
  return(list(m=m, Y=Y, S=S, X=X, subjects=ids))
}

## Function to build the data used in projections
build_pred_matrices<-function(df, future_years, subject_id, y_cols, s_cols,
                            fixed_vars, centers, stretches, levs,
                            half_life = 2.0){
  df <- as.data.frame(df)
  future_years <- sort(future_years)
  H <- length(future_years)
  
  ids <- sort(unique(df[[subject_id]]))
  n <- length(ids)
  last_year_raw <- max(df$year_raw)
  
  fixed_formula_text<-paste0("~ ", paste0(fixed_vars, collapse=" + "))
  fixed_formula<-as.formula(fixed_formula_text)
  
  X <- vector("list", n)
  Y <- vector("list", n)
  S <- vector("list", n)
  m <- NULL
  
  for (i in seq_along(ids)) {
    id <- ids[i]
    subj <- df[df[[subject_id]] == id, ]
    subj <- subj[order(subj$year_raw), ]
    last_row <- subj[which.max(subj$year_raw), , drop = FALSE]
    
    # recency-weighted estimate of this pitcher's future `stuff`, computed on
    # the already-scaled column so it's directly comparable/insertable
    w <- 0.5 ^ ((last_year_raw - subj$year_raw) / half_life)
    stuff_future_scaled <- sum(w * subj[[s_cols]]) / sum(w)
    
    # un-scale age back to real units so we can advance it by exactly the
    # number of years out, then re-scale with the SAME centers/stretches
    # used at training time
    last_age_raw <- last_row$age * stretches["age"] + centers["age"]
    
    panel <- last_row[rep(1, 1 + H), , drop = FALSE]
    rownames(panel) <- NULL
    panel$year_raw <- c(last_row$year_raw, future_years)
    
    age_raw_seq <- last_age_raw + (panel$year_raw - last_row$year_raw)
    panel$age  <- (age_raw_seq - centers["age"]) / stretches["age"]
    panel$year <- (panel$year_raw - centers["year"]) / stretches["year"]
    
    panel[[s_cols]] <- c(last_row[[s_cols]], rep(stuff_future_scaled, H))
    panel[[subject_id]] <- factor(id, levels = levs[[subject_id]])
    
    # future outcome counts are, by definition, unobserved
    panel[y_cols] <- NA_real_
    panel[1, y_cols] <- last_row[y_cols]
    panel$m <- c(last_row$m, rep(NA_real_, H))
    
    X[[i]] <- model.matrix(fixed_formula, data = panel)
    S[[i]] <- as.matrix(panel[, s_cols])
    Y[[i]] <- as.matrix(panel[, y_cols])
    m <- rbind(m, panel$m)
  }
  
  list(m = m, Y = Y, S = S, X = X, subjects = ids, future_years = future_years)
}

# Model fitting functions ------------------------------------------------------
## Function to find the MLE
MN_MLE<-function(n, m, X, Y, S, rho, B, lam=1, 
                 autoregress=TRUE, exogenous=TRUE,
                 eps=1e-3, max_iters=1000){
  k_minus_1<-ncol(Y[[1]])
  d<-1
  r<-1
  rho_new<-rho
  B_new<-B
  
  pi_hat<-calc_pi(n, m, Y, S, X, rho, B, k_minus_1, autoregress=autoregress,
                  exogenous=exogenous)
  
  ll<-c(sum(sapply(1:n, function(i)
    sum(
      Y[[i]] * log(pmax(pi_hat[[i]], 1e-12)),
      na.rm = TRUE
    ) +
      sum(m[i,] - rowSums(Y[[i]])) *
      log(pmax(1 - rowSums(pi_hat[[i]]), 1e-12))
  )))
  
  while(d>eps){
    ## Calculate autoregressive stuff slope
    if(autoregress){
      if(exogenous){
        grad<-nabla_rho(n, S, Y, m, pi_hat)
        hess<-nabla_2_rho(n, S, Y, m, pi_hat)
        
        ## Perform update step
        rho_new<-rho-solve(hess)%*%grad
      } else{
        grad<-nabla_rho_matrix(n, S, Y, m, pi_hat)
        hess<-nabla_2_rho_matrix(n, S, Y, m, pi_hat)
        
        ## Perform update step
        rho_new<-rho-matrix(solve(hess)%*%grad, ncol=k_minus_1, nrow=k_minus_1)
        
        pi_r<-calc_pi(n, m, Y, S, X, rho_new, B, k_minus_1,
                        autoregress=autoregress, exogenous=exogenous)
        
        ll_r<-sum(sapply(1:n, function(i)
          sum(Y[[i]]*log(pmax(pi_r[[i]], 1e-12)), na.rm=TRUE) +
            sum(m[i,]-rowSums(Y[[i]]))*log(pmax(1-rowSums(pi_r[[i]]), 1e-12))
        ))
        
        ll<-c(ll, ll_r)
        
      }
      
    }
    
    ## Calculate exogenous coefficients
    if(exogenous){
      grad<-nabla_B(n, X, Y, m, B, pi_hat)
      hess<-nabla_2_B(n, X, Y, m, B, pi_hat)
      
      ## Perform update step
      penalty<-diag(lam, ncol=length(B), nrow=length(B))
      step<-solve(hess-penalty)%*%(grad-lam*B)
      
      step_scale<-1
      repeat{
        B_try<-B-step_scale*step
        pi_try<-calc_pi(n, m, Y, S, X, rho_new, B_try, k_minus_1,
                        autoregress=autoregress, exogenous=exogenous)
        ll_try<-sum(sapply(1:n, function(i)
          sum(Y[[i]]*log(pmax(pi_try[[i]], 1e-12)), na.rm=TRUE) +
            sum(m[i,]-rowSums(Y[[i]]))*log(pmax(1-rowSums(pi_try[[i]]), 1e-12))
        ))
        if(ll_try>=ll[r] || step_scale<1e-6) break
        step_scale<-step_scale/2
      }
      
      ## Finalize step-halving components
      B_new<-B_try
      pi_hat<-pi_try
      ll<-c(ll, ll_try)
    }
    
    
    ## Determine new error
    d<-sqrt(sum((B_new - B)^2)+sum((rho-rho_new)^2))
    
    ## Update parameters and iteration count
    B<-B_new
    rho<-rho_new
    r<-r+1
    
    ## Stop if running too long
    if(r>=max_iters){
      break
    }
  }
  cat("========= MN MLE FIT =========\n")
  cat("Initial Log Likelihood: ", ll[1], "\n")
  cat("  Final log Likelihood: ", ll[length(ll)], "\n")
  cat("      Total iterations: ", r, "\n")
  cat("      Convergence norm: ", d, "\n")
  cat("  Attained convergence: ", d<eps, "\n")
  return(list(B=B, rho=rho, loglikelihood=ll, iters=r, conv_criterion=d))
}

## Function to evaluate model diagnostics, useful for cross-validation of lambda
mod_diags<-function(n, m, X, Y, S, rho, B, 
                    autoregress=TRUE, exogenous=TRUE, eval_last_only=FALSE){
  K<-ncol(Y[[1]])
  
  pi_hat<-calc_pi(n, m, Y, S, X, rho, B, K, autoregress=autoregress,
                  exogenous=exogenous)
  
  if(eval_last_only){
    J <- nrow(Y[[1]])
    observed <- do.call(rbind, lapply(1:n, function(i) Y[[i]][J,] / m[i,J]))
    predicted <- do.call(rbind, lapply(1:n, function(i) pi_hat[[i]][J,]))
  } else {
    observed <- do.call(rbind, lapply(1:n, function(i)
      Y[[i]] / m[i, ]
    ))
    
    predicted <- do.call(rbind, lapply(1:n, function(i)
      pi_hat[[i]]
    ))
  }
  
  ## Overall observed vs predicted rates
  obs_k_rate<-mean(observed[,1])
  pred_k_rate<-mean(predicted[,1])
  obs_bb_rate<-mean(observed[,2])
  pred_bb_rate<-mean(predicted[,2])
  
  cat("\n==== MODEL FIT ====\n")
  cat("Observed K rate:  ", obs_k_rate, "\n")
  cat("Predicted K rate: ", pred_k_rate, "\n")
  cat("Observed BB rate:  ", obs_bb_rate, "\n")
  cat("Predicted BB rate: ", pred_bb_rate, "\n")
  
  ## RMSE
  k_RMSE<-sqrt(mean((observed[,1] - predicted[,1])^2))
  bb_RMSE<-sqrt(mean((observed[,2] - predicted[,2])^2))
  cat("\nK RMSE:  ", k_RMSE, "\n")
  cat("BB RMSE: ", bb_RMSE, "\n")
  
  ## Correlation between observed and predicted rates
  k_corr<-cor(observed[,1], predicted[,1])
  bb_corr<-cor(observed[,2], predicted[,2])
  cat("\nK correlation:  ", k_corr, "\n")
  cat("BB correlation: ", bb_corr, "\n")
  
  ## Simple calibration plot
  p1<-ggplot()+
    geom_point(aes(x=observed[,1], y=predicted[,1]))+
    geom_abline(intercept=0, slope=1)+
    labs(title="Strikeout Calibration")+xlab("Observed K Rate")+ylab("Predicted K Rate")
  
  p2<-ggplot()+
    geom_point(aes(x=observed[,2], y=predicted[,2]))+
    geom_abline(intercept=0, slope=1)+
    labs(title="Walk Calibration")+xlab("Observed BB Rate")+ylab("Predicted BB Rate")
  
  return(list(obs_k_rate=obs_k_rate, pred_k_rate=pred_k_rate, 
              obs_bb_rate=obs_bb_rate, pred_bb_rate=pred_bb_rate,
              k_RMSE=k_RMSE, bb_RMSE=bb_RMSE, 
              k_corr=k_corr, bb_corr=bb_corr,
              k_plot=p1, bb_plot=p2))
}

## Function to perform cross-validation to select shrinkage parameter
MN_CV<-function(df, subject_id, y_cols, s_cols, fixed_vars, 
                lambdas=c(0.001, 0.01, 0.1, 1, 10, 20, 50, 75, 100),
                autoregress=TRUE, exogenous=TRUE, max_iters=1000){
  cv_results <- data.frame(
    lambda = lambdas,
    k_RMSE = NA_real_,
    bb_RMSE = NA_real_,
    k_corr = NA_real_,
    bb_corr = NA_real_
  )
  
  train_df <- df %>% filter(round(year_raw) <= 2023)
  
  test_ids <- intersect(
    unique(df$pitcher_id[round(df$year_raw) == 2023]),
    unique(df$pitcher_id[round(df$year_raw) == 2024])
  )
  
  test_df <- df %>%
    filter(pitcher_id %in% test_ids, round(year_raw) %in% c(2023, 2024)) %>%
    arrange(pitcher_id, year_raw)
  
  # Build training and validation matrices
  train_data <- build_matrices(
    train_df, subject_id, y_cols, s_cols, fixed_vars
  )
  
  K<-ncol(train_data$Y[[1]])
  
  test_data <- build_matrices(
    test_df, subject_id, y_cols, s_cols, fixed_vars
  )
  
  n_train <- length(train_data$subjects)
  n_test  <- length(test_data$subjects)
  
  train_pos <- setNames(seq_along(train_data$subjects), train_data$subjects)
  
  for(l in seq_along(lambdas)){
    
    cat("\n====================================\n")
    cat("Lambda:", lambdas[l], "\n")
    cat("====================================\n")
    
    # Initialize coefficient estimates
    rho_init<-rep(1, K)
    B_init <- rep(0, ncol(train_data$X[[1]]) * K)
    
    # Fit using 2020-2023
    fit <- MN_MLE(
      n = n_train,
      m = train_data$m,
      X = train_data$X,
      Y = train_data$Y,
      S = train_data$S,
      rho = rho_init,
      B = B_init,
      lam = lambdas[l],
      autoregress = autoregress,
      exogenous = exogenous,
      max_iters = max_iters
    )
    
    rho_test <- fit$rho
    
    # Evaluate ONLY on the 2024 row of each test panel — row 1 (2023) exists
    # solely to supply the lag-1 term.
    # (renamed from `diag` to `cv_diag` so it can't shadow base::diag())
    cv_diag <- mod_diags(
      n = n_test,
      m = test_data$m,
      X = test_data$X,
      Y = test_data$Y,
      S = test_data$S,
      rho = rho_test,
      B = fit$B,
      autoregress = autoregress,
      exogenous = exogenous,
      eval_last_only = TRUE
    )
    
    cv_results$k_RMSE[l]<-cv_diag$k_RMSE
    cv_results$bb_RMSE[l]<-cv_diag$bb_RMSE
    cv_results$k_corr[l]<-cv_diag$k_corr
    cv_results$bb_corr[l]<-cv_diag$bb_corr
  }
  
  # Print results
  print(cv_results)
  
  # Optimal lambdas
  best_k <- cv_results$lambda[which.min(cv_results$k_RMSE)]
  best_bb <- cv_results$lambda[which.min(cv_results$bb_RMSE)]
  
  cat("\n====================================\n")
  cat("Optimal lambda for K RMSE: ", best_k, "\n")
  cat("Optimal lambda for BB RMSE:", best_bb, "\n")
  cat("====================================\n")
  return(cv_results)
}

