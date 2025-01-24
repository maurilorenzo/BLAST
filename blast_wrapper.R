library(abind)
library(LaplacesDemon)
library(MASS)
library(matrixStats)

library(Rcpp)
library(RcppArmadillo)
sourceCpp('helper_functions.cpp')


pre_estimate_factors <- function(Y, k=NA, q_s=NA, k_max=50, tau=0.1, flag_svd_y=1, flag_svd_p=1,
                                 svd_cpp=F) {
  
  ns <- sapply(Y, function(x) nrow(x))
  S <- length(Y); p <- ncol(Y[[1]])
  Vs_tilde <- list()
  Vs_tilde_outer <- array(dim=c(p, p, S))
  
  qs <- rep(0, S)
  if(! is.na(k)){
    qs <- k + q_s
  }
  
  for(s in 1:S){
    n <- nrow(Y[[s]])
    if(is.na(k)){
      est_k <- estimate_latent_dimension(Y[[s]], k_max)
      qs[s] <- est_k$k_hat
      s_Y_s <- est_k$svd_Y
    } else{
      s_Y_s <- svd(Y[[s]], nu=k_max, nv=k_max)
    }
    Vs_tilde <-  s_Y_s$v[, 1:qs[s]]
    Vs_tilde_outer[,,s] <- Vs_tilde %*% t(Vs_tilde)
  }
  
  V_tilde_outer_mean <- apply(Vs_tilde_outer, c(1,2), mean)
  
  if(svd_cpp){
    s_L <- compute_svd_cpp(V_tilde_outer_mean, flag=flag_svd_p)
  } else {
    s_L <- svd(V_tilde_outer_mean, nu=k_max, nv=k_max)
  }
  if(is.na(k)){
    k <- sum((s_L$d>1-tau))
    q_s <- qs - k
    q_s[q_s<1] <- 1
  }
  V_bar <-  s_L$u[, 1:k]
  V_bar_outer <- V_bar %*% t(V_bar)
  V_bar_outer_perp <- diag(1, p, p) - V_bar_outer
  
  Y_shared <- matrix(NA, nrow=0, ncol=p)
  Y_hats_perp <- list()
  Fs <- list()
  Ms_2 <- list()
  
  
  for(s in 1:S){
    n <- nrow(Y[[s]])
    Y_hats_perp[[s]] <- Y[[s]] %*%  V_bar_outer_perp
    Fs[[s]] <- svd(Y_hats_perp[[s]])$u[,1:q_s[s]] * sqrt(n)
    QFs <- diag(1, n , n)  - Fs[[s]] %*% t(Fs[[s]]) / n 
    Ys_shared <- QFs %*%  Y[[s]]
    Y_shared <- rbind(Y_shared, Ys_shared)
    Ms_2[[s]] <- svd(Ys_shared)$u[,1:k] * sqrt(n)
  }
  if(svd_cpp){
    s_Y_shared <- compute_svd_cpp(Y_shared, flag=flag_svd_y)
  } else {
    s_Y_shared <- svd(Y_shared, nu=k, nv=k)
  }
  M_joint <- s_Y_shared$u[,1:k]*sqrt(sum(ns))
  prec_joint <- s_Y_shared$v[,1:k] %*% diag(s_Y_shared$d[1:k]) / sqrt(sum(ns))
  return(list(M=M_joint, Fs=Fs, prec=prec_joint, k=k, q_s=q_s, Y_shared=Y_shared,
              Ms_2 = Ms_2, V_bar_outer_perp=V_bar_outer_perp, P_tilde=V_tilde_outer_mean))
}

fit_blast <- function(
    Y, k = NA, q_s = NA, k_max = 30, n_MC = 100, cc_Lambda = T, cc_Gamma = T,  sample_outer_product = T, 
    subsample_index = 1:100, svd_cpp = T, svd_econ = T, tau = 0.1){
  
  ### INPUT ###
  # Y - list of matrices - list of length S where each element of is the a data matrix (n_s * p)
  # k - int - number of shared latent factors (if NA it's estimated)
  # q_s - vector(S) of int - number of study-specific latent factors
  # n_MC - int - number of Monte Carlo samples
  # cc_Lambda - bool (default to TRUE) - whether to apply coverage correction to Lambda 
  # cc_Gamma - bool (default to TRUE) - whether to apply coverage correction to the Gamma_s's 
  # sample_outer_product - bool (default to TRUE) - whether to sample outer products of loading matrices
  # subsample_index - vector of int - subset of indexes for which the LR components of the covariance are sampled (used only if sample_outer_product = TRUE) 
  # svd_cpp  bool (default to TRUE) - whether to use cpp implementation of SVD
  # svd_econ  bool (default to TRUE) - whether to use econ cpp implementation of SVD when p>1000 (used if svd_cpp = TRUE)
  # tau - double (default to 0.1) - hyperpameter used in shared number of factors selection (see Section 2.4 of the paper)
  
  ### OUTPUT ###
  # Lambda_outer_mean - matrix(p, p) - posterior mean of outer product of Lambda
  # Lambda_outer_mean - array(p, p, S) - posterior means of outer product of the Gamma_s's
  # Sigma_2s_samples - matrix(p, n_MC) - posterior samples of residual error variances
  # rho_Lambda - double - coverage correction factor for Lambda
  # rho_Gammas - vec(p) - coverage correction factors for the Gamma_s's
  # Ms - list(S) of Matrices - estimates of study-shared latent factors
  # Fs - list(S) of Matrices - estimates of study-specific latent factors
  # Lambda_mean - matrix(p, k) - point estimate for Lambda 
  # Gammas_mean - array(p, max(q_s), S) - point estimates for the Gamma_s's
  # P_tilde - matrix(p, p) - matrix defined in Section 2.1 of the main paper
  # Lambda_samples - array(p, k, n_MC) - posterior samples for Lambda
  # Gammas_samples - array(p, max(q_s), n_MC, S) - posterior samples for the Gamma_s's
  # Lambda_outer_samples - array(length(subsample_index), length(subsample_index), n_MC) - posterior samples for outer product of Lambda (returned only if sample_outer_product=TRUE)
  # Gammas_outer_samples - array(length(subsample_index), length(subsample_index), n_MC, S) - posterior samples for outer product of the Gamma_s's (returned only if sample_outer_product=TRUE)
  
  
  S <- length(Y); p <- ncol(Y[[1]])
  ns <- sapply(Y, function(x) nrow(x))
  n_joint <- sum(ns)
  p_sample <- length(subsample_index)
  flag_svd_y=1; flag_svd_p=1
  if(svd_econ & p > 1000){
    flag_svd_y=3; flag_svd_p=4
  }
  
  factors_estimates <- pre_estimate_factors(
    Y, k=k, q_s=q_s, k_max=k_max, flag_svd_y=flag_svd_y, flag_svd_p=flag_svd_p,
    svd_cpp=svd_cpp, tau=tau
    )
  k <- factors_estimates$k
  q_s <- factors_estimates$q_s
  
  Sigma_2s_samples <- array(0, dim=c(p, n_MC, S))
  Gammas_samples <- array(0, dim=c(p, max(q_s), n_MC, S))
  Gammas_outer_mean <- array(0, dim=c(p, p, S))
  if(sample_outer_product){
    Gammas_outer_samples <- array(0, dim=c(p_sample, p_sample, n_MC, S))
    
  }
  Gammas_mean <- array(0, dim=c(p, max(q_s), S))
  rho_Gammas <- rep(1, S)
  print(paste0('number of shared latent factors: ', k))
  print(paste0('number of study-specific latent factors: ', q_s))

  M_joint <- factors_estimates$M
  Fs <- factors_estimates$Fs
  Ms <- compute_Ms(M_joint, ns, S)
  Y_joint <-  do.call(abind, c(Y, along = 1))
  Y_perp_joint <- Y_joint - 1/n_joint* tcrossprod(M_joint) %*% Y_joint
  
  Y_shared <- factors_estimates$Y_shared
  Y_joint <-  Y_shared
  
  L_js_s <- compute_L_js(Y_joint, M_joint); V_js<- compute_V_js(Y_joint, M_joint)
  tau_2_hat_s <- 1/(k) * sum(L_js_s) / sum(V_js)
  print('fitting shared loadings')
  mu_Lambda <- t(Y_joint) %*% M_joint / (n_joint + 1/tau_2_hat_s)
  mu_Lambda_samples <- array(rep(mu_Lambda, n_MC), dim = c(n_MC, p, k))
  rho_Lambda <- 1

  
  if(cc_Lambda){
    B <- compute_B_Lambda(mu_Lambda, V_js)
    rho_Lambda <- mean(B[lower.tri(B,diag=TRUE)])
    print(rho_Lambda)
  }
  
  Lambda_samples <- sample_Lambda(Y_joint, mu_Lambda, tau_2_hat_s, M_joint, 
                                  gamma_0=1, delta_0=1, N_mc=n_MC, rho=rho_Lambda,
                                  subsample_index=subsample_index, sample_outer_product=FALSE)
  Lambda_outer_mean <-  Lambda_samples$Lambda_outer_mean

  if(sample_outer_product){
    Lambda_outer_samples <- array(NA, dim=c(p_sample, p_sample, n_MC))
    Lambda_outer_samples[,,] <- sample_Lambda_outer(Lambda_samples$Lambda_samples[subsample_index,,])
  }
  
  index_start <- 1
  for(s in 1:S){
    print(paste0('fitting study-specific loadings for study ', s))
    n_s <- ns[s]
    index_finish <- index_start + ns[s]-1
    Y_perp_s <- Y_perp_joint[index_start:index_finish, ]
    index_start <- index_finish+1
    
    if(q_s[s]>0){
      Y_perp <- Y[[s]] - Ms[[s]] %*% t(mu_Lambda)
      L_js_s <- compute_L_js(Y_perp, Fs[[s]])
      V_js_s <- compute_V_js(Y_perp, Fs[[s]])
      tau_2_hat_s <- 1/(q_s[s]) * sum(L_js_s)/sum(V_js_s)
      print(tau_2_hat_s)
      n_s <- nrow(Y_perp)
      mu_js_s <- t(Y_perp) %*% Fs[[s]] /(n_s +  1/tau_2_hat_s)
      
      if(cc_Gamma){
        B <- compute_B_Gamma(mu_Lambda, mu_js_s, V_js_s)
        rho_Gammas[s] <- mean(B[lower.tri(B,diag=TRUE)])
        print(rho_Gammas[s])
      }
      Gamma_s_fable <- sample_Gamma(
        Y[[s]], tau_2_hat_s, Fs[[s]], Lambda_samples$Lambda_samples, Ms[[s]],
        Lambda_samples$Sigmas_samples, gamma_0=1, delta_0=1, N_mc=n_MC, rho=rho_Gammas[s],
        subsample_index=subsample_index-1)
      Gammas_mean[,1:q_s[s],s] <- mu_js_s
      Ms_2 <- factors_estimates$Ms_2[[s]]
      y_perp <- Y[[s]] - 1/n_s* tcrossprod(Ms_2) %*% Y[[s]]
      F_ <- svd(Y_perp_s)$u[,1:q_s[s]] *sqrt(n_s)
      Gammas_samples[,1:q_s[s],,s] <- Gamma_s_fable$Gamma_samples[,,]
      Gammas_outer_mean[,,s] <- tcrossprod(mu_js_s) + rho_Gammas[s]^2 * diag(as.vector(V_js))/n_s
      if(sample_outer_product){
        Gammas_outer_samples[,,,s] <- Gamma_s_fable$Gamma_outer_samples
      }
    }
  }
  
  
  output <- list(Lambda_outer_mean = Lambda_outer_mean,
                 Gammas_outer_mean = Gammas_outer_mean,
                 Sigma_2s_samples = Lambda_samples$Sigmas_samples,
                 rho_Lambda=rho_Lambda, rho_Gammas=rho_Gammas,
                 Ms = Ms, Fs = Fs,
                 Lambda_mean = mu_Lambda,
                 Gammas_mean = Gammas_mean, 
                 P_tilde = factors_estimates$P_tilde,
                 Lambda_samples = Lambda_samples$Lambda_samples,
                 Gammas_samples = Gammas_samples)
  if(! sample_outer_product){
    return(output)
  }
  
  output$Lambda_outer_samples = Lambda_outer_samples
  output$Gammas_outer_samples = Gammas_outer_samples
  
  return(output)
}


compute_Ms_cross <- function(M, n_s, S){
  index_start <- 1
  index_finish <- n_s[1]
  Ms_cross <- list()
  Ms_cross[[1]] <- crossprod(M[1:n_s[1],], M[1:n_s[1],])
  for(s in 2:S){
    index_start <- index_finish + 1
    index_finish <- index_finish + n_s[s]
    Ms_cross[[s]] <- crossprod(M[index_start:index_finish,], 
                               M[index_start:index_finish,])
    
  }
  return(Ms_cross)
}

compute_Ms <- function(M, n_s, S){
  index_start <- 1
  index_finish <- n_s[1]
  Ms <- list()
  Ms[[1]] <-M[1:n_s[1],]
  for(s in 2:S){
    index_start <- index_finish + 1
    index_finish <- index_finish + n_s[s]
    Ms[[s]] <- M[index_start:index_finish,]
    
  }
  return(Ms)
}

compute_tau_2 <- function(Y, M){
  p <- ncol(Y); k <- ncol(M)
  L_js_s <- compute_L_js(Y, M)
  V_js_s <- compute_V_js(Y, M)
  tau_2_hat <- 1/(k*p) * sum(L_js_s/V_js_s)
  return(tau_2_hat)
}

compute_jic <- function(Y, svd_Y, k){
  n <- nrow(Y); p <- ncol(Y) 
  minint <- min(n ,p)
  maxint <- max(n, p)
  M <- sqrt(n)*as.matrix(svd_Y$u[,1:k])
  Lambda <- 1/sqrt(n)* as.matrix(svd_Y$v[,1:k]) %*% diag(svd_Y$d[1:k], k, k)
  
  Y_hat <- tcrossprod(M, Lambda)
  sigma_sq_hat <- colMeans((Y - Y_hat)^2) # p * 1
  tausq_est <- (mean(colSums((Y_hat)^2) / sigma_sq_hat)) / (n * k);
  
  Lambda_est <- (sqrt(n) / (n + 1/tausq_est)) * svd_Y$v[,1:k] %*% diag(svd_Y$d[1:k], k, k);
  M_Lambda_est <- sqrt(n)* svd_Y$u[,1:k] %*% t(Lambda_est)
  res <- colSums((Y - M_Lambda_est)^2) 
  sigma_sq_hat <- res / n
  
  term1 <- -n*sum(log(sqrt(sigma_sq_hat))) - 0.5*sum(res/sigma_sq_hat)
  term1 <- term1 / (n*p)
  #print(k)
  jic <- -2*term1 + (k * maxint * log(minint) / (n*p))
  #print(jic)
  return(jic)
}

estimate_latent_dimension <- function(Y, k_max){
  svd_Y <- svd(Y)
  n <- nrow(Y)
  jics <- sapply(1:k_max, function(x) (compute_jic(Y, svd_Y, x)))
  plot(1:k_max, jics, type='l', xlab='k', ylab='jic', main='')
  print(paste('k_hat = ', which.min(jics)))
  return(list(k_hat = which.min(jics), jics=jics, svd_Y = svd_Y))
}

compute_L_js <- function(Y, U){
  R <- solve(t(U) %*% U) %*% t(U) %*% Y 
  return(colSums(R^2))
}

compute_V_js <- function(Y, U){
  R <- Y - U%*%solve(t(U) %*% U) %*%t(U) %*% Y
  return(colMeans(R^2))
}

