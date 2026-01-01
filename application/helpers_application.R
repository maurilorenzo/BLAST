# R SCRIPT WITH HELPER FUNCTIONS FOR APPLICATION CONTAINED IN SECTION 5 OF THE PAPER

if (!require("genefilter", quietly = TRUE)){
  BiocManager::install("genefilter")
}
if(! require(VIMSFA, quietly = TRUE)){
  install_github("blhansen/VI-MSFA")
}
library(Biobase)
library(circlize)
library(emdbook)
library(expm)
library(genefilter)
library(ggplot2)
library(mvtnorm)
library(pracma)
library(readxl)
library(sparsepca)
library(VIMSFA)


source('blast_wrapper.R')
source('helpers_competitors.R')


predict_factors <- function(Y_old, Lambda_old, Psi_old, coverage){
  k <- ncol(Lambda_old)
  Var_fact <- solve(diag(1, k, k) + t(Lambda_old) %*% diag(1/Psi_old) %*% Lambda_old)
  mean_fact <- Y_old %*% diag(1/Psi_old) %*% Lambda_old  %*% Var_fact 
  if(coverage){  eta <-mean_fact + mvtnorm::rmvnorm(nrow(Y_old), sigma=Var_fact)}
  output <- list(Etas_mean=mean_fact)
  if(coverage){output$Etas_sample=eta}
  return(output)
}

predict_y <- function(Y_old, Lambda_new, Lambda_old, Psi_new, Psi_old, coverage=F){
  Eta <- predict_factors(Y_old, Lambda_old, Psi_old, coverage)
  if(coverage){  Y_new <- Eta$Etas_sample %*% t(Lambda_new) +  mvtnorm::rmvnorm(nrow(Y_old), sigma=diag(Psi_new))}
  Y_new_mean <- Eta$Etas_mean %*% t(Lambda_new)
  output <- list(Y_mean=Y_new_mean)
  if(coverage){output$Y_sample=Y_new}
   return(output)
}

compute_oos_accuracy_blast <- function(blast_fit, Y_test, impute_factor_index_gene, coverage=F){
  n_MC <- dim(blast_fit$Lambda_samples)[3]
  n_test <- sapply(Y_test, nrow)
  p_test <- ncol(Y_test[[1]]) - length(impute_factor_index_gene)
  Y_preds <- list()
  Y_means <- list()
  #q_s <- sapply(blast_fit$Gammas_samples, function(x) dim(x)[2])
  for(s in 1:3){
    if(coverage){ Y_preds[[s]] <- array(0, dim=c(n_test[s], p_test, n_MC))}
    Y_means[[s]] <- array(0, dim=c(n_test[s], p_test))
    for(i in 1:n_MC){
      Loadings_est_gene <- cbind(blast_fit$Lambda_samples[,,i], 
                                 blast_fit$Gammas_samples[,,i,s])
      predict_i <- predict_y(
        Y_test_gene[[s]][,impute_factor_index_gene], Loadings_est_gene[-impute_factor_index_gene,], 
        Loadings_est_gene[impute_factor_index_gene,], 
        as.vector(blast_fit$Sigma_2s_samples[i,-impute_factor_index_gene]), 
        as.vector(blast_fit$Sigma_2s_samples[i, impute_factor_index_gene]),
        coverage=coverage)
      if(coverage){Y_preds[[s]][,,i] <- predict_i$Y_sample}
      Y_means[[s]][,] <- Y_means[[s]][,] + predict_i$Y_mean
    } 
    Y_means[[s]][,] <- Y_means[[s]][,] / n_MC
  }
  
  blast_gene_coverage <- rep(NA, 3)
  blast_gene_mses <- array(NA, dim=c(3, p_test))
  blast_gene_mses_norm <- array(NA, dim=c(3, p_test))
  
  
  for(s in 1:3){
    print(s)
    if(coverage){
      Y_qs <- apply(Y_preds[[s]], c(1,2), function(x)(quantile(x, probs=c(0.025, 0.975))))
      blast_gene_coverage[s] <- mean(Y_qs[1,,]<Y_test_gene[[s]][,-impute_factor_index_gene] &
                                       Y_qs[2,,]>Y_test_gene[[s]][,-impute_factor_index_gene])
      print(blast_gene_coverage[s])
      
    }
    n_s <- nrow(Y_test_gene[[s]])
    blast_gene_mses[s,] <- colMeans( (Y_means[[s]] - Y_test_gene[[s]][,-impute_factor_index_gene])^2)
    blast_gene_mses_norm[s,] <- blast_gene_mses[s,]/ (colSds((Y_test_gene[[s]][,-impute_factor_index_gene])))^2 *  n_s / (n_s - 1)
    print( c( mean(blast_gene_mses_norm[s,]),
              median(blast_gene_mses_norm[s,]), 
              quantile(blast_gene_mses_norm[s,], probs=c(0.25, 0.75))
    ))
  }
  return(list(coverage=blast_gene_coverage, 
              mses = blast_gene_mses,
              mses_norm = blast_gene_mses_norm,
              Y_means=Y_means,
              Y_preds=Y_preds))
}

compute_samples_svi_Lambda_outer <- function(Lambda_mean, Lambda_vars, n_MC=100, 
                                             subsample_index=1:100, outer=T){
  p <- length(subsample_index); k <- ncol(Lambda_mean)
  Lambda_vars_root <- lapply(Lambda_vars, function(x) expm::sqrtm(x + 0.0000001*diag(1, k,k)))
  if(outer){
    samples_save <- array(NA, dim=c(n_MC, p, p))
  }
  else{
    samples_save <- array(NA, dim=c(n_MC, p, k))
  }
  for(i in 1:n_MC){
    Lambda_sample <- Lambda_mean[subsample_index, ]
    for(j in 1:p){
      j_p <- subsample_index[j]
      e_j <- rnorm(k) 
      Lambda_sample[j,] = Lambda_sample[j_p,] + t(Lambda_vars_root[[j_p]] %*% e_j)
    }
    if(outer){
      samples_save[i,,]  <- Lambda_sample %*% t(Lambda_sample)
    }
    else{
      samples_save[i,,]  <- Lambda_sample
    }
  }
  return(samples_save)
}

compute_samples_svi_Gammas_outer <- function(Gammas_means, Gammas_vars, n_MC=100,
                                             subsample_index=1:100, outer=T){
  S <- length(Gammas_means)
  p <- length(subsample_index)
  qs <- sapply(Gammas_means, ncol)
  if(outer){
    samples_save <- array(NA, dim=c(S, n_MC, p, p))
  } else {
    samples_save <- array(NA, dim=c(S, n_MC, p, max(qs)))
  }
  for(s in 1:S){
    samples_save[s,,,1:qs[s]] <- compute_samples_svi_Lambda_outer(
      Gammas_means[[s]], Gammas_vars[[s]], n_MC=n_MC, subsample_index=subsample_index,
      outer=outer)
  }
  return(samples_save) 
}


compute_oos_accuracy_vi <- function(vi_fit, Y_test, impute_factor_index_gene, n_MC=100, coverage=F){
  n_test <- sapply(Y_test, nrow)
  p_test <- ncol(Y_test[[1]]) - length(impute_factor_index_gene)
  p <- ncol(Y_test[[1]])
  qs <- sapply(vi_fit$mean_lambda_s, ncol)
  Y_preds <- list()
  Y_means <- list()
  Lambda_samples <- compute_samples_svi_Lambda_outer(
    vi_fit$mean_phi, vi_fit$var_phi, n_MC=n_MC, subsample_index=1:p, outer=F)
  Gammas_samples <- compute_samples_svi_Gammas_outer(
    vi_fit$mean_lambda_s, vi_fit$var_lambda_s, n_MC=n_MC, subsample_index=1:p, outer=F)
  for(s in 1:3){
    if(coverage){
      Y_preds[[s]] <- array(0, dim=c(n_test[s], p_test, n_MC))
    }
    Y_means[[s]] <- array(0, dim=c(n_test[s], p_test))
    
    for(i in 1:n_MC){
      Lambda_i <- Lambda_samples[i,,]
      Gamma_s_i <- Gammas_samples[s,i,,1:qs[s]]
      Loadings_i <- cbind(Lambda_i, Gamma_s_i)
      vars_i <- 1/sapply(vi_fit$rate_psi_s[[s]], 
                         function(x)(rgamma(1, shape=vi_fit$shape_psi_s[[s]], rate=x)))
      predict_i <- predict_y(
        Y_test_gene[[s]][,impute_factor_index_gene], Loadings_i[-impute_factor_index_gene,], 
        Loadings_i[impute_factor_index_gene,], vars_i[-impute_factor_index_gene], 
        vars_i[impute_factor_index_gene], coverage=coverage)
      if(coverage){Y_preds[[s]][,,i] <- predict_i$Y_sample}
      Y_means[[s]][,] <- Y_means[[s]][,] + predict_i$Y_mean
    } 
    Y_means[[s]][,] <- Y_means[[s]][,] / n_MC
  }
  
  blast_gene_coverage <- rep(NA, 3)
  blast_gene_mses <- array(NA, dim=c(3, p_test))
  blast_gene_mses_norm <- array(NA, dim=c(3, p_test))
  
  for(s in 1:3){
    print(s)
    if(coverage){
      Y_qs <- apply(Y_preds[[s]], c(1,2), function(x)(quantile(x, probs=c(0.025, 0.975))))
      blast_gene_coverage[s] <- mean(Y_qs[1,,]<Y_test_gene[[s]][,-impute_factor_index_gene] &
                                       Y_qs[2,,]>Y_test_gene[[s]][,-impute_factor_index_gene])
      print(blast_gene_coverage[s])
    }
    
    n_s <- nrow(Y_test_gene[[s]])
    blast_gene_mses[s,] <- colMeans( (Y_means[[s]] - Y_test_gene[[s]][,-impute_factor_index_gene])^2)
    blast_gene_mses_norm[s,] <- blast_gene_mses[s,] /  (colSds((Y_test_gene[[s]][,-impute_factor_index_gene])))^2 *  n_s / (n_s - 1)
    print( c( mean(blast_gene_mses_norm[s,]),
              median(blast_gene_mses_norm[s,]), 
              quantile(blast_gene_mses_norm[s,], probs=c(0.25, 0.75))
    ))
  }
  return(list(coverage=blast_gene_coverage, 
              mses = blast_gene_mses,
              mses_norm = blast_gene_mses_norm,
              Y_means = Y_means,
              Y_preds = Y_preds))
}


compute_oos_accuracy_ando_bai <- function(ando_bai_est, Y_test, impute_factor_index_gene, n_MC=100, coverage=F){
  n_test <- sapply(Y_test, nrow)
  p_test <- ncol(Y_test[[1]]) - length(impute_factor_index_gene)
  p <- ncol(Y_test[[1]])
  qs <- sapply(ando_bai_est$Gammas, ncol)
  Y_preds <- list()
  Y_means <- list()
  
  for(s in 1:3){
    Y_means[[s]] <- array(0, dim=c(n_test[s], p_test))
    Lambda_i <- ando_bai_est$Lambda_c
    Gamma_s  <- ando_bai_est$Gammas[[s]] #do.call(cbind, lapply(V_s, function(Vg) Vg[j, 1:k]))
    vars_i <- ando_bai_est$Sigmas[[s]]
    Loadings_i <- cbind(Lambda_i, Gamma_s)
    predict_i <- predict_y(
      Y_test_gene[[s]][,impute_factor_index_gene], Loadings_i[-impute_factor_index_gene,], 
      Loadings_i[impute_factor_index_gene,], vars_i[-impute_factor_index_gene], 
      vars_i[impute_factor_index_gene], coverage=F)
    Y_means[[s]][,] <- predict_i$Y_mean 
    if(coverage){
      Y_preds[[s]] <- array(0, dim=c(n_test[s], p_test, n_MC))
      for(i in 1:n_MC){
        predict_i <- predict_y(
          Y_test_gene[[s]][,impute_factor_index_gene], Loadings_i[-impute_factor_index_gene,], 
          Loadings_i[impute_factor_index_gene,], vars_i[-impute_factor_index_gene], 
          vars_i[impute_factor_index_gene], coverage=T)
        Y_preds[[s]][,,i] <- predict_i$Y_sample
      }
    }
    
  }
  
  blast_gene_coverage <- rep(NA, 3)
  blast_gene_mses <- array(NA, dim=c(3, p_test))
  blast_gene_mses_norm <- array(NA, dim=c(3, p_test))
  
  for(s in 1:3){
    print(s)
    if(coverage){
      Y_qs <- apply(Y_preds[[s]], c(1,2), function(x)(quantile(x, probs=c(0.025, 0.975))))
      blast_gene_coverage[s] <- mean(Y_qs[1,,]<Y_test_gene[[s]][,-impute_factor_index_gene] &
                                       Y_qs[2,,]>Y_test_gene[[s]][,-impute_factor_index_gene])
      print(blast_gene_coverage[s])
    }
    
    n_s <- nrow(Y_test_gene[[s]])
    blast_gene_mses[s,] <- colMeans( (Y_means[[s]] - Y_test_gene[[s]][,-impute_factor_index_gene])^2)
    blast_gene_mses_norm[s,] <- blast_gene_mses[s,] /  (colSds((Y_test_gene[[s]][,-impute_factor_index_gene])))^2 *  n_s / (n_s - 1)
    print( c( mean(blast_gene_mses_norm[s,]),
              median(blast_gene_mses_norm[s,]), 
              quantile(blast_gene_mses_norm[s,], probs=c(0.25, 0.75))
    ))
  }
  return(list(coverage=blast_gene_coverage, 
              mses = blast_gene_mses,
              mses_norm = blast_gene_mses_norm,
              Y_means = Y_means,
              Y_preds = Y_preds))
}



estimate_latent_dimensions <-  function(Y, tau=0.1){
  k_hats <- list()
  S <-length(Y)
  p <- ncol(Y[[1]])
  P_tilde <- matrix(0, p, p)
  k_max <- 0
  k_s <- c()
  for(s in 1:S){
    k_hats[[s]] <- estimate_latent_dimension(Y[[s]], k_max=100)
    V <-  k_hats[[s]]$svd_Y$v[,1: k_hats[[s]]$k_hat]; 
    P <- tcrossprod(V)
    P_tilde <- P_tilde + P
    k_max <- k_max + k_hats[[s]]$k_hat
    k_s[s] <- k_hats[[s]]$k_hat
  }
  P_tilde <- P_tilde / S
  s_tilde <- svd(P_tilde, nu=k_max, nv=k_max)
  k_0_hat <- max(which(s_tilde$d[1:(k_max)] > 1 - tau))
  return(list(k_0=k_0_hat, q_s = k_s - k_0_hat))
}


