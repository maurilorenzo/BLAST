
library(devtools)
if(! require(MSFA)){
  install_github("rdevito/MSFA")
}
if(! require(VIMSFA)){
  install_github("blhansen/VI-MSFA")
}


library(dplyr)
library(expm)
library(infinitefactor)
library(MSFA)
library(pracma)
library(readxl)
library(sparsepca)
library(VIMSFA)

source('blast_wrapper.R')
source('helpers_competitors.R')


################################################################################
# Helper functions for simulations setup
################################################################################

generate_data <- function(p=500, n_s=rep(500, 5), S=5, q_s=rep(5,5), k_0=10, seed=123, var='het'){
  set.seed(seed)
  sigma_loadings <- 0.5 # std of loadings
  pi_lambda = 0.5 # add sparsity (pi_lambda = prob of being non zero)
  l_sigma=0.5; u_sigma=5 # upper and lower bounds of residual variances
  Loadings <- matrix(rnorm(p*(k_0+sum(q_s)), 0, sigma_loadings), ncol=(k_0+sum(q_s)))
  Loadings <- Loadings * matrix(rbinom(p*(k_0+sum(q_s)), 1, pi_lambda), nrow=p)
  Lambda_0 <- Loadings[,1:k_0]
  Lambda_0_outer <- Lambda_0 %*% t(Lambda_0)
  
  Gammas_0 <- list()
  Gammas_0_outer <- list()
  Sigmas_0 <- list()
  LRs_0<-  list()
  Thetas_0 <- list()
  
  k_init <- k_0 +1
  
  for(s in 1:S){
    print(s)
    n <- n_s[s]
    Gamma_1 <- Loadings[, k_init:(k_init + q_s[s]-1)]
    k_init <- k_init + q_s[s]
    Gammas_0[[s]] <- Gamma_1
    Gammas_0_outer[[s]] <-  Gamma_1 %*% t(Gamma_1)
    LRs_0[[s]] <- Lambda_0_outer + Gammas_0_outer[[s]] 
    Sigmas_0[[s]] <- diag(runif(p, l_sigma, u_sigma), ncol=p, nrow=p)
    Thetas_0[[s]] <- LRs_0[[s]] + Sigmas_0[[s]]
  }
  
  
  ## generate data
  Etas_0 <- list()
  Phis_0 <- list()
  Y <- list()
  #set.seed(123)
  for(s in 1:S){
    n <- n_s[s]
    Etas_0[[s]] <- matrix(rnorm(n*k_0), ncol=k_0)
    Phis_0[[s]] <- matrix(rnorm(n*q_s[s]), ncol=q_s[s])
    Y[[s]] <- Etas_0[[s]] %*% t(Lambda_0)+  Phis_0[[s]] %*% t(Gammas_0[[s]])
    if(var == 'het'){
      Y[[s]] <- Y[[s]] + mvrnorm(n, rep(0,p), Sigmas_0[[s]])
    }
    else{
      Y[[s]] <- Y[[s]] + mvrnorm(n, rep(0,p), Sigmas_0[[1]])
    }
  }
  return(list(Y=Y, Lambda_0_outer=Lambda_0_outer, Gammas_0_outer=Gammas_0_outer,
              Etas_0=Etas_0, Phis_0=Phis_0))
}



print_metrics <- function(df_results, S=5){
  df <- as.matrix(df_results)
  n_tries <- nrow(df_results)
  print(c('error', 'coverage', 'length'))
  print('Lambda')
  print(colMeans(df[,2:4]))
  print(colSds(df[,2:4])/sqrt(n_tries))
  print('Gammas')
  print(c(mean(df[,5:(4+S)]), mean(df[,(5+S):(4+2*S)]),  mean(df[,(5+2*S):(4+3*S)])) )
  print(c(sd(df[,5:(4+S)]), sd(df[,(5+S):(4+2*S)]),  sd(df[,(5+2*S):(4+3*S)]))/sqrt(n_tries))
  print('Etas & Phis')
  print(c(mean(df[,(5+3*S):(4+4*S)]), mean(df[,(5+4*S):(4+5*S)])))
  print(c(sd(df[,(5+3*S):(4+4*S)]), sd(df[,(5+4*S):(4+5*S)]))/sqrt(n_tries))
  print('time')
  print(mean(df[, ncol(df)]))
  print(sd(df[, ncol(df)])/sqrt(n_tries))
}


################################################################################
# Helper functions for BLAST
################################################################################

compute_metrics_blast <- function(blast_fit, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0,
                                  n_MC=500, subsample_index=1:100){
  Lambda_0_outer_sub <- Lambda_0_outer[subsample_index, subsample_index]
  output_blast <- c()
  output_blast[1] = norm(tcrossprod(blast_fit$Lambda_mean) - Lambda_0_outer, type='F') / norm(Lambda_0_outer, type='F')
  pr_hom_blast_estimate_Lambda_qs <- apply(blast_fit$Lambda_outer_samples, c(1,2), 
                                           function(x)(quantile(x, probs=c(0.025, 0.975))))
  output_blast[2] = mean((pr_hom_blast_estimate_Lambda_qs[1,,] < Lambda_0_outer_sub) &
                           (pr_hom_blast_estimate_Lambda_qs[2,,] > Lambda_0_outer_sub))
  output_blast[3] = mean(pr_hom_blast_estimate_Lambda_qs[2,,] - pr_hom_blast_estimate_Lambda_qs[1,,])
  for(s in 1:S){
    print(s)
    # reconstruction study specific LR part
    output_blast[s+3] = norm(tcrossprod(blast_fit$Gammas_mean[,,s]) - Gammas_0_outer[[s]], type='F') /
      norm(Gammas_0_outer[[s]], type='F')
  }
  
  pr_hom_blast_estimate_Lambda_qs <- apply(
    blast_fit$Lambda_outer_samples[, subsample_index, subsample_index ], c(2,3), 
    function(x) (quantile(x, probs=c(0.025, 0.975)))
  )
  for(idx in 1:S){
    pr_hom_blast_estimate_Gamma_qs <- apply(
      blast_fit$Gammas_outer_samples[,subsample_index,subsample_index,idx], c(2,3), 
      function(x) (quantile(x, probs=c(0.025, 0.975)))
    )
    output_blast[S+3+ idx] = mean((pr_hom_blast_estimate_Gamma_qs[1,,] < Gammas_0_outer[[idx]][subsample_index,subsample_index ])
                                  & (pr_hom_blast_estimate_Gamma_qs[2,,] > Gammas_0_outer[[idx]][subsample_index,subsample_index]))
    output_blast[2*S+3+ idx] = mean(pr_hom_blast_estimate_Gamma_qs[2,,] - pr_hom_blast_estimate_Gamma_qs[1,,] )
  }
  
  for(s in 1:S){
    mean_f_1 <- blast_fit$Ms[[s]]
    p_1 <- procrustes(mean_f_1, Etas_0[[s]])
    n_1 <- nrow(mean_f_1); k_1 <- ncol(mean_f_1)
    print(norm(mean_f_1 - Etas_0[[s]] %*% p_1$Q, type='F') / norm(Etas_0[[s]], type='F'))
    output_blast[3*S+3+ s] = norm(mean_f_1 - Etas_0[[s]] %*% p_1$Q, type='F') / sqrt(n_1*k_1)
  }
  
  for(s in 1:S){
    mean_f_1 <-  blast_fit$Fs[[s]]
    dim(mean_f_1)
    n_1 <- nrow(mean_f_1); q_1 <- ncol(mean_f_1)
    p_1 <- procrustes(mean_f_1, Phis_0[[s]])
    print(norm(mean_f_1 - Phis_0[[s]] %*% p_1$Q, type='F') / norm(Phis_0[[s]], type='F'))
    output_blast[4*S+3+ s] = norm(mean_f_1 - Phis_0[[s]] %*% p_1$Q, type='F') / sqrt(n_1*q_1)
    
  }
  return(output_blast)
}


################################################################################
# Helper functions for VI-MSFA
################################################################################

compute_samples_svi_Lambda_outer <- function(Lambda_mean, Lambda_vars, n_MC=100, subsample_index=1:100){
  p <- length(subsample_index); k <- ncol(Lambda_mean)
  Lambda_vars_root <- lapply(Lambda_vars, function(x) expm::sqrtm(x + 0.0000001*diag(1, k,k)))
  Lambda_outer_samples <- array(NA, dim=c(n_MC, p, p))
  for(i in 1:n_MC){
    Lambda_sample <- Lambda_mean[subsample_index, ]
    for(j in 1:p){
      j_p <- subsample_index[j]
      e_j <- rnorm(k) 
      Lambda_sample[j,] = Lambda_sample[j_p,] + t(Lambda_vars_root[[j_p]] %*% e_j)
    }
    Lambda_outer_samples[i,,]  <- Lambda_sample %*% t(Lambda_sample)
  }
  return(Lambda_outer_samples)
}

compute_samples_svi_Gammas_outer <- function(Gammas_means, Gammas_vars, n_MC=100, subsample_index=1:100){
  S <- length(Gammas_means)
  p <- length(subsample_index)
  Gammas_outer_samples <- array(NA, dim=c(S, n_MC, p, p))
  for(s in 1:S){
    Gammas_outer_samples[s,,,] <- compute_samples_svi_Lambda_outer(Gammas_means[[s]], Gammas_vars[[s]], n_MC=n_MC, subsample_index=subsample_index)
  }
  return(Gammas_outer_samples) 
}



compute_metrics_vi <- function(vi_fit, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0,
                               n_MC=500, subsample_index=1:100){
  vi_phi <- vi_fit$mean_phi
  vi_lambda_s <- vi_fit$mean_lambda_s
  vi_psi_s <- vi_fit$mean_psi_s
  vi_var_phi <- vi_fit$var_phi
  vi_var_lambda_s <- vi_fit$var_lambda_s
  output_vi <- c()
  output_vi[1] = norm(vi_phi %*% t(vi_phi) - Lambda_0_outer, type='F') / norm(Lambda_0_outer, type='F')
  #print('sample lambda')
  vi_Lambda_samples <- compute_samples_svi_Lambda_outer(vi_phi, vi_var_phi, n_MC=500, 
                                                        subsample_index=subsample_index)
  vi_Lambda_qs <- apply(vi_Lambda_samples, c(2,3), function(x)(quantile(x, probs=c(0.025, 0.975))))
  output_vi[2] = mean((Lambda_0_outer_sub > vi_Lambda_qs[1,,]) & (Lambda_0_outer_sub < vi_Lambda_qs[2,,]))
  output_vi[3] = mean(vi_Lambda_qs[2,,] - vi_Lambda_qs[1,,])
  rm(vi_Lambda_samples)
  rm(vi_Lambda_qs)
  
  for(i in 1:S){
    output_vi[i+3] = norm(vi_lambda_s[[i]] %*% t(vi_lambda_s[[i]]) - Gammas_0_outer[[i]], type='F') / norm(Gammas_0_outer[[i]], type='F')
  }
  
  #print('sample gammas')
  vi_Gammas_outer_samples <- compute_samples_svi_Gammas_outer(vi_fit$mean_lambda_s, vi_fit$var_lambda_s, n_MC=100, subsample_index=subsample_index)
  vi_Gammas_outer_qs <- apply(vi_Gammas_outer_samples, c(1, 3,4), function(x)(quantile(x, probs=c(0.025, 0.975))))
  for(s in 1:S){
    output_vi[S+3+s] = mean((Gammas_0_outer[[s]][subsample_index, subsample_index] > vi_Gammas_outer_qs[1,s,,]) & (Gammas_0_outer[[s]][subsample_index,subsample_index]  < vi_Gammas_outer_qs[2,s,,]))
    output_vi[2*S+3+s] =  mean(vi_Gammas_outer_qs[2,s,,] -  vi_Gammas_outer_qs[1,s,,])
  }
  
  for(s in 1:S){
    #print('procrustes error')
    mean_f_1 <- do.call(rbind, lapply(vi_fit$mean_f[[s]], as.vector))
    n_s <- nrow(mean_f_1)
    k_0 <-  ncol(mean_f_1)
    p_1 <- procrustes(mean_f_1, Etas_0[[s]])
    output_vi[3*S+3+s] = norm(mean_f_1 - Etas_0[[s]] %*% p_1$Q, type='F') / sqrt(n_s*k_0)
    mean_l_1 <- do.call(rbind, lapply(vi_fit$mean_l[[s]], as.vector))
    q_s <-  ncol(mean_l_1)
    p_1 <- procrustes(mean_l_1, Phis_0[[s]])
    output_vi[4*S+3+s] = norm(mean_l_1 - Phis_0[[s]] %*% p_1$Q, type='F') /sqrt(n_s*q_s)
  }
  
  rm(vi_Gammas_outer_samples)
  return(output_vi)
}

################################################################################
# Helper functions for MSFA and BMSFA
################################################################################
compute_outer_samples_msfa <- function(bmsfa_fit, subsample_index=1:p){
  n_MC <-  dim(bmsfa_fit$Phi)[3]
  p <- dim(bmsfa_fit$Phi)[1]
  p_sample <- length(subsample_index)
  S <- length(bmsfa_fit$Lambda)
  Lambda_outer_samples <- array(0, dim=c(n_MC, p_sample, p_sample))
  q_s <- sapply(bmsfa_fit$Lambda, function(x) (dim(x)[2]))
  Gammas_outer_samples <- array(0, dim=c(S, n_MC, p_sample, p_sample))
  
  for(it in 1:n_MC) {
    Lambda_outer_samples[it,,] <- tcrossprod(bmsfa_fit$Phi[subsample_index,,it])
    for(s in 1:S){
      Gammas_outer_samples[s,it,,] <- tcrossprod(bmsfa_fit$Lambda[[s]][subsample_index,,it]) 
    }
  }
  return(list(Lambda_outer_samples=Lambda_outer_samples, 
              Gammas_outer_samples=Gammas_outer_samples))
}


full_conditional_factors <- function(Y, Loadings, sigma){
  k <- ncol(Loadings)
  n <- nrow(Y)
  var_ <- solve(diag(1,k,k) + t(Loadings) %*% diag(1/sigma) %*% Loadings)
  mean_ <- Y %*% diag(1/sigma) %*% Loadings %*% var_
  sample <- mean_ + mvrnorm(n, mu=rep(0, k) ,Sigma=var_)
  return(sample)
}

compute_factor_samples_bmsfa <- function(bmsfa_fit, Y){
  p <- dim(bmsfa_fit$Phi)[1]
  k <-  dim(bmsfa_fit$Phi)[2]
  n_MC <- dim(bmsfa_fit$Phi)[3]
  
  S <- length(bmsfa_fit$Lambda)
  eta_samples <- list();   phi_samples <- list()
  q_s <- sapply(bmsfa_fit$Lambda, function(x) (dim(x)[2]))
  #Gammas_outer_samples <- array(0, dim=c(S, n_MC, p , p ))
  for(s in 1:S){
    n <- nrow(Y[[s]])
    eta_samples[[s]] <- array(NA, dim=c(n, k, n_MC))
    phi_samples[[s]] <- array(NA, dim=c(n, q_s[s], n_MC))
    
    for(it in 1:n_MC) {
    Lambda_sample <- bmsfa_fit$Phi[,,it]
    Gamma_sample <- bmsfa_fit$Lambda[[s]][,,it]
    sigma_sample <- as.vector(bmsfa_fit$psi[[s]][,,it])
    Loadings <- cbind(Lambda_sample, Gamma_sample)
    factors <- full_conditional_factors(Y[[s]], Loadings, sigma_sample)
    eta_samples[[s]][,,it] <- factors[,1:k]
    phi_samples[[s]][,,it] <- factors[,-c(1:k)]
    }
  }
  return(list(eta_samples=eta_samples, 
              phi_samples=phi_samples))
} 

compute_mean_factors_bmsfa <- function(bmsfa_fit){
  S <- length(bmsfa_fit$f_s)
  Etas <- list(); Phis <- list()
  for(s in 1:S){
    Etas_s <- matchalign(bmsfa_fit$f_s[[s]], bmsfa_fit$Phi)$eta
    Etas[[s]] <- Reduce("+", Etas_s) / length(Etas_s)
    Phis_s <- matchalign(bmsfa_fit$l_s[[s]], bmsfa_fit$Lambda[[s]])$eta
    Phis[[s]] <-  Reduce("+", Phis_s) / length(Phis_s)
  }
  
  return(list(Etas=Etas, Phis=Phis))
}

convert_array_to_list <- function(array) {
  n_MC <- dim(array)[3]
  return(lapply(seq_len(n_MC), function(i) array[, , i]))
}

matchalign <- function(factors, loadings){
  factors_list <- convert_array_to_list(factors)
  loadings_list <- convert_array_to_list(loadings)
  out <- jointRot(loadings_list, factors_list)
  return(out)
} 





compute_metrics_bmsfa <- function(bmsfa_fit, Lambda_0_outer, Gammas_0_outer, Etas_0, 
                                  Phis_0,  subsample_index=1:100){
  
  p <- dim(bmsfa_fit$Phi)[1]
  n_MC <- dim(bmsfa_fit$Phi)[3]
  Lambda_0_outer_sub <- Lambda_0_outer[subsample_index, subsample_index]
  bmsfa_samples <- compute_outer_samples_msfa(bmsfa_fit, subsample_index=1:p)
  bmsfa_Lambda_outer_mean <- apply(bmsfa_samples$Lambda_outer_samples, c(2, 3), mean)
  bmsfa_Gammas_outer_mean <- apply(bmsfa_samples$Gammas_outer_samples, c(1, 3, 4), mean)
  bmsfa_Lambda_outer_qs <- apply(bmsfa_samples$Lambda_outer_samples[,subsample_index, subsample_index], c(2, 3), function(x)
    (quantile(x, probs=c(0.025, 0.975))))
  bmsfa_Gammas_outer_qs <- apply(bmsfa_samples$Gammas_outer_samples[,,subsample_index, subsample_index], c(1, 3, 4), function(x)
    (quantile(x, probs=c(0.025, 0.975))))
  
  output_bmsfa <- c()
  output_bmsfa[1] = norm(bmsfa_Lambda_outer_mean - Lambda_0_outer, type='F') / norm(Lambda_0_outer, type='F')
  output_bmsfa[2] = mean((Lambda_0_outer_sub > bmsfa_Lambda_outer_qs[1,,]) & 
                          (Lambda_0_outer_sub < bmsfa_Lambda_outer_qs[2,,]))
  output_bmsfa[3] = mean(bmsfa_Lambda_outer_qs[2,,] - bmsfa_Lambda_outer_qs[1,,])
  rm(bmsfa_Lambda_outer_qs)

  for(i in 1:S){
    output_bmsfa[i+3] = norm(bmsfa_Gammas_outer_mean[i,,] - Gammas_0_outer[[i]], type='F') / norm(Gammas_0_outer[[i]], type='F')
  }
  for(s in 1:S){
    output_bmsfa[S+3+s] = mean((Gammas_0_outer[[s]][subsample_index, subsample_index] > bmsfa_Gammas_outer_qs[1,s,,]) &
                                 (Gammas_0_outer[[s]][subsample_index, subsample_index]  < bmsfa_Gammas_outer_qs[2,s,,]))
    output_bmsfa[2*S+3+s] =  mean(bmsfa_Gammas_outer_qs[2,s,,] -  bmsfa_Gammas_outer_qs[1,s,,])
  }

  rm(bmsfa_Gammas_outer_qs)
  rm(bmsfa_samples)
  
  factors_samples <- compute_factor_samples_bmsfa(bmsfa_fit, Y)
  bmsfa_fit$f_s <- factors_samples$eta_samples
  bmsfa_fit$l_s <- factors_samples$phi_samples
  aligned_factors <- compute_mean_factors_bmsfa(bmsfa_fit)
  mean_eta <- aligned_factors$Etas; mean_phi <- aligned_factors$Phis; 

  for(s in 1:S){
    mean_f_1 <- mean_eta[[s]]
    n_s <- nrow(mean_f_1)
    k_0 <-  ncol(mean_f_1)
    p_1 <- procrustes(mean_f_1, Etas_0[[s]])
    output_bmsfa[3*S+3+s] = norm(mean_f_1 - Etas_0[[s]] %*% p_1$Q, type='F') / sqrt(n_s*k_0)
    mean_l_1 <- mean_phi[[s]]
    q_s <-  ncol(mean_l_1)
    p_1 <- procrustes(mean_l_1, Phis_0[[s]])
    output_bmsfa[4*S+3+s] = norm(mean_l_1 - Phis_0[[s]] %*% p_1$Q, type='F') /sqrt(n_s*q_s)
  }
  return(output_bmsfa)
}

compute_Thompson_factor_scores <- function(Y, Lambda, Psi){
  k <- ncol(Lambda)
  Z <- Y %*% diag(1/Psi) %*% Lambda %*% solve(t(Lambda) %*% diag(1/Psi) %*% Lambda + diag(1, k,k))
  return(Z)
}


compute_mean_factors_msfa <- function(msfa_fit, Y){
  S <- length(msfa_fit$Lambda_s)
  k <- ncol(msfa_fit$Phi)
  etas <- list()
  phis <- list()
  for(s in 1:S){
    Loadings <- cbind(msfa_fit$Phi, msfa_fit$Lambda_s[[s]])
    sigma <- as.vector(msfa_fit$psi_s[[s]])
    factors <- compute_Thompson_factor_scores(Y[[s]], Loadings, sigma)
    etas[[s]] <- factors[,1:k]; 
    phis[[s]] <- factors[,-c(1:k)]; 
  }
  return(list(etas=etas, phis=phis))
}

compute_metrics_msfa <- function(msfa_est, Y, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0, 
                                 n_MC=n_MC, subsample_index=subsample_index){

  output_msfa <- c()
  output_msfa[1] = norm(tcrossprod(msfa_est$Phi) - Lambda_0_outer, type='F') / norm(Lambda_0_outer, type='F')
  output_msfa[2] = NA
  output_msfa[3] = NA

  for(s in 1:S){
    output_msfa[s+3] = norm(tcrossprod(msfa_est$Lambda_s[[s]]) - Gammas_0_outer[[s]], type='F') / norm(Gammas_0_outer[[s]], type='F')
  }
  for(s in 1:S){
    output_msfa[S+3+s] = NA    
    output_msfa[2*S+3+s] =  NA
  }
  
  mean_factors <- compute_mean_factors_msfa(msfa_fit, Y)
  mean_eta <- mean_factors$etas; mean_phi <- mean_factors$phis; 
  
  for(s in 1:S){
    Loadings_s <- cbind(msfa_est$Phi, msfa_est$Lambda_s[[s]])
    k_0 <- ncol(msfa_est$Phi); q_s <- ncol(msfa_est$Lambda_s[[s]])
    factors_s <- compute_Thompson_factor_scores(Y[[s]], Loadings_s, as.vector(msfa_est$psi_s[[s]]))
    mean_f_1 <- factors_s[,1:k_0]
    n_s <- nrow(mean_f_1)
    p_1 <- procrustes(mean_f_1, Etas_0[[s]])
    output_msfa[3*S+3+s] = norm(mean_f_1 - Etas_0[[s]] %*% p_1$Q, type='F') / sqrt(n_s*k_0)
    mean_l_1 <- factors_s[,-c(1:k_0)]
    q_s <-  ncol(mean_l_1)
    p_1 <- procrustes(mean_l_1, Phis_0[[s]])
    output_msfa[4*S+3+s] = norm(mean_l_1 - Phis_0[[s]] %*% p_1$Q, type='F') /sqrt(n_s*q_s)
  }
  return(output_msfa) 
}

################################################################################
# Helper functions for spectral estimator by Bai & Ng
################################################################################

compute_metrics_spectral <- function(spectral_fit, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0,
                                  subsample_index=1:100){
  Lambda_0_outer_sub <- Lambda_0_outer[subsample_index, subsample_index]
  output_spectral <- c()
  S <- length(Gammas_0_outer)
  output_spectral[1] = norm(tcrossprod(spectral_fit$Lambda_c) - Lambda_0_outer, type='F') / norm(Lambda_0_outer, type='F')
  for(s in 1:S){
    print(s)
    # reconstruction study specific LR part
    output_spectral[s+3] = norm(tcrossprod(spectral_fit$Gammas[[s]]) - Gammas_0_outer[[s]], type='F') /
      norm(Gammas_0_outer[[s]], type='F')
  }
  
  
  for(s in 1:S){
    mean_f_1 <- spectral_fit$Ms[[s]]
    p_1 <- procrustes(mean_f_1, Etas_0[[s]])
    n_1 <- nrow(mean_f_1); k_1 <- ncol(mean_f_1)
    print(norm(mean_f_1 - Etas_0[[s]] %*% p_1$Q, type='F') / norm(Etas_0[[s]], type='F'))
    output_spectral[3*S+3+ s] = norm(mean_f_1 - Etas_0[[s]] %*% p_1$Q, type='F') / sqrt(n_1*k_1)
  }
  
  for(s in 1:S){
    mean_f_1 <-  spectral_fit$Fs[[s]]
    dim(mean_f_1)
    n_1 <- nrow(mean_f_1); q_1 <- ncol(mean_f_1)
    p_1 <- procrustes(mean_f_1, Phis_0[[s]])
    print(norm(mean_f_1 - Phis_0[[s]] %*% p_1$Q, type='F') / norm(Phis_0[[s]], type='F'))
    output_spectral[4*S+3+ s] = norm(mean_f_1 - Phis_0[[s]] %*% p_1$Q, type='F') / sqrt(n_1*q_1)
    
  }
  return(output_spectral)
}

