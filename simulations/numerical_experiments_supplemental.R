### R SCRIPT TO REPLICATE THE NUMERICAL EXPERIMENTS IN THE SUPPLEMENTARY MATERIAL

source('simulations/helpers_numerical_experiments.R')


S <- 3; k_0 <- 4; q_s <- rep(3, S)

p <- 150; n_s <- rep(200, S); var <- 'hom'

test_svi <- T; test_cavi <- T; test_blast <- F; test_msfa <- F; test_bmsfa <- F

n_sim <- 50
sim <- 1

df_svi <- data.frame()
df_cavi <- data.frame()
df_blast <- data.frame()
df_bmsfa <- data.frame()
df_msfa <- data.frame()


for(sim in 1:n_sim){
  print(sim)
  data_sim <- generate_data(p=p, n_s=n_s, S=S, q_s=q_s, k_0=k_0, seed=sim, var=var)
  Y=data_sim$Y; Lambda_0_outer=data_sim$Lambda_0_outer; Gammas_0_outer=data_sim$Gammas_0_outer;
  Etas_0=data_sim$Etas_0; Phis_0=data_sim$Phis_0
  subsample_index <- 1:min(p, 100)
  Lambda_0_outer_sub <- Lambda_0_outer[subsample_index, subsample_index]
  n_MC = 500
  if(test_svi){
    set.seed(123)
    ptm <- proc.time()
    svi_est <- svi_msfa(Y, K=k_0, J_s=q_s, batch_prop = 0.2, max_iter = 1000, scale=F)
    svi_time =  proc.time() - ptm
    output_svi <- compute_metrics_vi(svi_est, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0,
                                     n_MC=n_MC, subsample_index=subsample_index)
    rm(svi_est)
    output_svi[5*S +4] = svi_time[3]
    df_svi <- rbind(df_svi, c(1, output_svi))
  }
  if(test_cavi){
    set.seed(123)
    ptm <- proc.time()
    cavi_est <- cavi_msfa(Y, k_0, q_s, scale=F, max_iter=1000)
    cavi_time <-  proc.time() - ptm
    output_cavi <- compute_metrics_vi(cavi_est, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0, n_MC=n_MC,
                                      subsample_index=subsample_index)
    rm(cavi_est)
    output_cavi[5*S +4] = cavi_time[3]
    df_cavi <- rbind(df_cavi, c(2, output_cavi))
  }
  if(test_bmsfa){
    set.seed(123)
    ptm <- proc.time()
    bmsfa_est <- sp_msfa(Y, k=k_0, j_s=q_s, trace=FALSE, control=list(nrun=10000, burn=5000, thin=10))
    bmsfa_time <-  proc.time() - ptm
    output_bmsfa <- compute_metrics_bmsfa(bmsfa_est, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0, 
                                          subsample_index=subsample_index)
    rm(bmsfa_est)
    output_bmsfa[5*S +4] = bmsfa_time[3]
    df_bmsfa <- rbind(df_bmsfa, c(4, output_bmsfa))
  }
  if(test_msfa){
    set.seed(123)
    ptm <- proc.time()
    msfa_start_value <- start_msfa(Y, k=k_0, j_s=q_s)
    msfa_est <-  ecm_msfa(Y, msfa_start_value, trace=FALSE, corr=FALSE)
    msfa_time <-  proc.time() - ptm
    output_msfa <- compute_metrics_msfa(msfa_est, Y, Lambda_0_outer, Gammas_0_outer, Etas_0, 
                                        Phis_0, subsample_index=subsample_index)
    rm(msfa_est)
    output_msfa[5*S +4] = msfa_time[3]
    df_msfa <- rbind(df_msfa, c(5, output_msfa))
  }
  if(test_blast){
    set.seed(123)
    ptm <- proc.time()
    blast_est <- fit_blast(Y, k=NA, q_s=NA, n_MC=n_MC, k_max=20,
                           subsample_index=subsample_index,
                           svd_cpp=T, svd_econ=T, tau=0.2) 
    blast_time =  proc.time() - ptm
    output_blast <- compute_metrics_blast(blast_est, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0,
                                          n_MC=n_MC, subsample_index=subsample_index)
    output_blast[5*S + 4] = blast_time[3]
    df_blast <- rbind(df_blast, c(3, output_blast))
  }
  
}


names <- c('method', 'L_fr', 'L_cov', 'L_len', paste0("G", 1:S, "_fr"), paste0("G", 1:S, "_cov"), 
           paste0("G", 1:S, "_len"), paste0("Eta", 1:S, "_fr"), paste0("Phi", 1:S, "_fr"),
           'time')
names(df_svi) <- names
names(df_cavi) <- names
names(df_bmsfa) <- names
names(df_msfa) <- names
names(df_blast) <- names

print_metrics(df_svi, S=3)
print_metrics(df_cavi, S=3)
print_metrics(df_bmsfa, S=3)
print_metrics(df_msfa, S=3)
print_metrics(df_blast, S=3)
