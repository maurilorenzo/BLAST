### R SCRIPT TO REPLICATE THE NUMERICAL EXPERIMENTS IN SECTION 4

source('simulations/helpers_numerical_experiments.R')


S <- 5 # number of studies
k_0 <- 10 # latent dimension of shared component
q_s <- rep(5, S) # latent dimension of study-specific component

# change scenario from 1 to 8 for all configurations considered
scenario <- 1
p <- 5000 # outcome dimension 
n_s <- rep(1000, S)# sample sizes
var <- 'hom'

if(scenario > 5){
  var <- 'het'
}
if((scenario %% 2) == 0){
  n_s <- rep(1000, S)
}
if(scenario %in% c(3,4,7,8)){
  p <- 5000
}

df_svi <- data.frame()
df_cavi <- data.frame()
df_blast <- data.frame()

test_svi <- T; test_cavi <- T; test_blast <- T
n_sim <- 5
for(sim in 4:n_sim){
  print(sim)
  data_sim <- generate_data(p=p, n_s=n_s, S=S, q_s=q_s, k_0=k_0, seed=sim, var=var)
  Y=data_sim$Y; Lambda_0_outer=data_sim$Lambda_0_outer; Gammas_0_outer=data_sim$Gammas_0_outer;
  Etas_0=data_sim$Etas_0; Phis_0=data_sim$Phis_0
  subsample_index <- 1:100
  Lambda_0_outer_sub <- Lambda_0_outer[subsample_index, subsample_index]
  
  if(test_svi){
    ptm <- proc.time()
    svi_est <- svi_msfa(Y, K=k_0, J_s=q_s, batch_prop = 0.2, max_iter = 1000, scale=F)
    svi_time =  proc.time() - ptm
    print(svi_time[3])
    output_svi <- compute_metrics_vi(svi_est, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0,
                                     n_MC=50, subsample_index=subsample_index)
    rm(svi_est)
    output_svi[5*S +4] = svi_time[3]
    df_svi <- rbind(df_svi, c(1, output_svi))
  }
  if(test_cavi){
    ptm <- proc.time()
    cavi_est <- cavi_msfa(Y, k_0, q_s, scale=F, max_iter=1000)
    cavi_time <-  proc.time() - ptm
    print(cavi_time[3])
    output_cavi <- compute_metrics_vi(cavi_est, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0, n_MC=50,
                                      subsample_index=subsample_index)
    rm(cavi_est)
    output_cavi[5*S +4] = cavi_time[3]
    df_cavi <- rbind(df_cavi, c(2, output_cavi))
  }
  if(test_blast){
    ptm <- proc.time()
    blast_est <- fit_blast(Y, k=NA, q_s=NA, n_MC=500, k_max=20, subsample_index=subsample_index,
                           svd_cpp=T, svd_econ=T) 
    blast_time =  proc.time() - ptm
    blast_time[3]
    output_blast <- compute_metrics_blast(blast_est, Lambda_0_outer, Gammas_0_outer, Etas_0, Phis_0,
                                          n_MC=500, subsample_index=subsample_index)
    output_blast[5*S + 4] = blast_time[3]
    df_blast <- rbind(df_blast, c(3, output_blast))
  }
}


names <- c('method', 'L_fr', 'L_cov', 'L_len', paste0("G", 1:S, "_fr"), paste0("G", 1:S, "_cov"), 
           paste0("G", 1:S, "_len"), paste0("Eta", 1:S, "_fr"), paste0("Phi", 1:S, "_fr"),
           'time')
names(df_svi) <- names
names(df_cavi) <- names
names(df_blast) <- names

print_metrics(df_svi[,])
print_metrics(df_cavi)
print_metrics(df_blast)

