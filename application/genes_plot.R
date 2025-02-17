# R SCRIPT TO MAKE PLOTS FOR THE IN APPLICATION IN SECTION 5 OF THE PAPER, 
# WHICH ARE REPORTED IN THE SUPPLEMENTAL MATERIAL


source('application/helpers_application.R')
cutoff <- 0.75
source('application/gene_preprocessing.R')

set.seed(123)
n_MC<- 100
ptm <- proc.time()
blast_gene_analysis <- fit_blast(
  Y, k=k_0_hat, q_s=q_s_hat, n_MC=n_MC,
  sample_outer_product = T, k_max=120, subsample_index = 1:1000) 
blast_gene_analysis_time_1 <- (proc.time() - ptm)[3]
blast_gene_analysis_time_1

subsample_plot <- 1:1000
Lambda_outer_blast_gene_sel <- blast_gene_analysis$Lambda_outer_mean[subsample_plot, subsample_plot]
Gammas_outer_blast_gene_sel <- blast_gene_analysis$Gammas_outer_mean[subsample_plot, subsample_plot,]
Lambda_outer_blast_gene_sel_qs <- apply(blast_gene_analysis$Lambda_outer_samples,
                                        c(1,2),
                                        function(x)(quantile(x, probs=c(0.025, 0.975))))
Lambda_outer_blast_gene_sel_est <- Lambda_outer_blast_gene_sel
Lambda_outer_blast_gene_sel_est[Lambda_outer_blast_gene_sel_qs[1,,]<0 &
                                  Lambda_outer_blast_gene_sel_qs[2,,]>0] = 0
C2=cov2cor(Lambda_outer_blast_gene_sel); 
d <- as.dist(1-abs(C2)) 
set.seed(123)
hc2=cluster::agnes(d)

Gammas_outer_blast_gene_sel_qs <- apply(blast_gene_analysis$Gammas_outer_samples,
                                        c(1, 2, 4),
                                        function(x)(quantile(x, probs=c(0.025, 0.975))))
Gammas_outer_blast_gene_sel_est <- Gammas_outer_blast_gene_sel
S <- 3
for(s in 1:S){
  print(mean(Gammas_outer_blast_gene_sel_qs[1,,,s]<0 & 
               Gammas_outer_blast_gene_sel_qs[2,,,s]>0))
  Gammas_outer_blast_gene_sel_est[,,s][Gammas_outer_blast_gene_sel_qs[1,,,s]<0 & 
                                         Gammas_outer_blast_gene_sel_qs[2,,,s]>0] = 0
  print(mean(Gammas_outer_blast_gene_sel_est[,,s] == 0))
}
mean(Gammas_outer_blast_gene_sel_est==0)

################################################################################
# FIG 1 supp
################################################################################
library(ggplot2)
custom_palette <- colorRampPalette(c("#AA4499", "white", "#117733"))
blast_gene_sigmas_est <- colMeans(blast_gene_analysis$Sigma_2s_samples)

# reconstructed correlation
corr_1 <- cov2cor(Gammas_outer_blast_gene_sel[hc2$order,hc2$order, 1] + 
                    Lambda_outer_blast_gene_sel[hc2$order,hc2$order] + 
                    diag(blast_gene_sigmas_est[hc2$order])
                  )
Loadings_1 <- blast_gene_analysis$Lambda_outer_samples+ blast_gene_analysis$Gammas_outer_samples[,,,1]
Loadings_outer_blast_gene_sel_qs <- apply(Loadings_1,
                                        c(1,2),
                                        function(x)(quantile(x, probs=c(0.025, 0.975))))
corr_1[Loadings_outer_blast_gene_sel_qs[1,hc2$order,hc2$order]<0 &
         Loadings_outer_blast_gene_sel_qs[2,hc2$order,hc2$order]>0] = 0

corr_data <- as.data.frame(
  as.table(corr_1))
corr_data$study = 'GSE15907'

for(s in 2:S){
  corr_s <- cov2cor(Gammas_outer_blast_gene_sel[hc2$order,hc2$order, s] + 
                      Lambda_outer_blast_gene_sel_est[hc2$order,hc2$order] + 
                      diag(blast_gene_sigmas_est[hc2$order]))
  Loadings_1 <- blast_gene_analysis$Lambda_outer_samples+ blast_gene_analysis$Gammas_outer_samples[,,,s]
  Loadings_outer_blast_gene_sel_qs <- apply(Loadings_1,
                                            c(1,2),
                                            function(x)(quantile(x, probs=c(0.025, 0.975))))
  
  corr_s[Loadings_outer_blast_gene_sel_qs[1,hc2$order,hc2$order]<0 &
           Loadings_outer_blast_gene_sel_qs[2,hc2$order,hc2$order]>0] = 0
  corr_data_s <- as.data.frame(
    as.table(corr_s))
  if(s==2){
    #corr_data_s$study = paste0('Study ',s)
    corr_data_s$study = 'GSE37448'
  }
  else{
    corr_data_s$study = 'GSE109125'
  }
  corr_data <- rbind(corr_data, corr_data_s) 
}
corr_shared <- cov2cor(Lambda_outer_blast_gene_sel[hc2$order,hc2$order])
corr_shared[Lambda_outer_blast_gene_sel_qs[1,hc2$order,hc2$order]<0 &
               Lambda_outer_blast_gene_sel_qs[2,hc2$order,hc2$order]>0] = 0
corr_data_s <- as.data.frame(as.table(corr_shared))
#corr_data_s$study = 'Shared'
#corr_data <- rbind(corr_data, corr_data_s) 
lim_min = -1; lim_max= 1
corr_plots_2 <- ggplot(corr_data, aes(Var1, Var2, fill=Freq)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = custom_palette(100), 
    limits = c(lim_min, lim_max),
    name = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_text(size = 14))+
  labs(x="", y="", fill="") +
  facet_wrap(~ study, nrow=1)
corr_plots_2
png(file = "fig/corr_plots_1000.png", width = 1000, height = 300) 
corr_plots_2
dev.off()
dev.off()

# empirical correlation
emp_corr_1 <- cor(Y[[1]][, subsample_plot][, hc2$order])
emp_corr_data <- as.data.frame(as.table(emp_corr_1))
emp_corr_data$study = 'GSE15907'
for(s in 2:S){
  emp_corr_s <- cor(Y[[s]][, subsample_plot][, hc2$order])
  emp_corr_data_s <- as.data.frame(
    as.table(emp_corr_s))
  if(s==2){
    #corr_data_s$study = paste0('Study ',s)
    emp_corr_data_s$study = 'GSE37448'
  }
  else{
    emp_corr_data_s$study = 'GSE109125'
  }
  #emp_corr_data_s$study = paste0('Study ',s)
  emp_corr_data <- rbind(emp_corr_data, emp_corr_data_s) 
}
lim_min = -1; lim_max= 1

custom_palette <- colorRampPalette(c("#AA4499", "white", "#117733"))
emp_corr_plot<- ggplot(emp_corr_data, aes(Var1, Var2, fill=Freq)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = custom_palette(100),
    limits = c(lim_min, lim_max),
    name = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_text(size = 14))+
  labs(x="", y="", fill="") +
  facet_wrap(~ study, nrow=1)
#dev.new()
emp_corr_plot

png(file = "fig/emp_corr_plot_1000.png", width = 1000, height = 300) 
emp_corr_plot
dev.off()

################################################################################
# FIG 2 supp
################################################################################

# shared vs study-specifi LR
corr_1 <- cov2cor(Gammas_outer_blast_gene_sel[hc2$order,hc2$order, 1])
corr_1[Gammas_outer_blast_gene_sel_qs[1,hc2$order, hc2$order,1]<0 & 
  Gammas_outer_blast_gene_sel_qs[2, hc2$order, hc2$order,1]>0] = 0
corr_data <- as.data.frame(as.table(corr_1))
corr_data$study = 'GSE15907'
for(s in 2:S){
  corr_s <- cov2cor(Gammas_outer_blast_gene_sel[hc2$order,hc2$order, s])
  corr_s[Gammas_outer_blast_gene_sel_qs[1,hc2$order,hc2$order,s]<0 & 
           Gammas_outer_blast_gene_sel_qs[2,hc2$order,hc2$order,s]>0] = 0
  corr_data_s <- as.data.frame(
    as.table(corr_s))
  if(s==2){
    corr_data_s$study = 'GSE37448'
  }
  else{
    corr_data_s$study = 'GSE109125'
  }
  corr_data <- rbind(corr_data, corr_data_s) 
}
corr_shared <- cov2cor(Lambda_outer_blast_gene_sel[hc2$order,hc2$order])
corr_shared[Lambda_outer_blast_gene_sel_qs[1,hc2$order,hc2$order]<0 &
              Lambda_outer_blast_gene_sel_qs[2,hc2$order, hc2$order]>0] = 0
corr_data_s <- as.data.frame(as.table(corr_shared))
corr_data_s$study = 'Shared'
corr_data <- rbind(corr_data, corr_data_s) 
lim_min = -1; lim_max= 1
custom_palette <- colorRampPalette(c("#AA4499", "white", "#117733"))
corr_plots_3 <- ggplot(corr_data, aes(Var1, Var2, fill=Freq)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = custom_palette(100), 
    limits = c(lim_min, lim_max), 
    name = ""
  ) +
  theme_minimal() + 
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    strip.text = element_text(size = 14))+
  labs(x="", y="", fill="") +
  facet_wrap(~ study, nrow=1)

corr_plots_3
png(file = "fig/corr_shared_vs_study_plot_1000.png", width = 1300, height = 300) 
corr_plots_3
dev.off()


################################################################################
# FIG 3 supp
################################################################################

sigma_2s_blast_gene <- colMeans(blast_gene_analysis$Sigma_2s_samples)
shared_cov <- Lambda_outer_blast_gene_sel + diag(sigma_2s_blast_gene[subsample_plot])
cor_mean <- cov2cor(shared_cov)
cor_mean[Lambda_outer_blast_gene_sel_qs[1,,]<0 & Lambda_outer_blast_gene_sel_qs[2,,]>0] = 0

rownames(cor_mean)=colnames(cor_mean)=genes.use[subsample_plot]
cor_mean[abs(cor_mean)<.5]=0
mean(cor_mean==0)
indices.connected.genes=which(rowSums( cor_mean!=0)-1 > 10)
length(indices.connected.genes)
subgenes=genes.use[indices.connected.genes]

C2=cor_mean[subgenes,subgenes]; 
d2 <- as.dist(1-abs(C2)) 
hc2=cluster::agnes(d2)
indices.connected.genes.ordered=hc2$order
subgenes.connected.ordered=subgenes[indices.connected.genes.ordered]



# this creates csv files for GEPHI (Bastian et al. 2009)

cov_graph <- Lambda_outer_blast_gene_sel[indices.connected.genes, indices.connected.genes]
cor_graph <- cov2cor(cov_graph)
cor_graph[Lambda_outer_blast_gene_sel_qs[1,indices.connected.genes,indices.connected.genes]<0 & 
           Lambda_outer_blast_gene_sel_qs[2,indices.connected.genes,indices.connected.genes]>0] = 0
nodes <- subgenes.connected.ordered
edges <- data.frame()
threshold <- 0.5

for (i in 1:(nrow(cor_graph)-1)) {
  for (j in (i+1):ncol(cor_graph)) {
    if (cor_graph[i,j]  > threshold) {
      new_edge <- data.frame(Source = nodes[i], Target = nodes[j], Weight = cor_graph[i,j], Type='Undirected')
      edges <- rbind(edges, new_edge)
    }
  }
}

nodes_df <- data.frame(list(Id=seq(1, 209), Label=subgenes.connected.ordered))
write.csv(edges, 'application/fig/edges_df.csv')
write.csv(nodes_df, 'application/fig/nodes_df.csv')



# adapting code from https://github.com/noirritchandra/SUFA/blob/main/vignettes/Genedata_application.R
est_partcor_plot_subgenes_connected_ordered=C2[subgenes.connected.ordered,subgenes.connected.ordered] 
lim_min <- min(est_partcor_plot_subgenes_connected_ordered)
lim_max <- max(est_partcor_plot_subgenes_connected_ordered)
c(lim_min, lim_max)
col_viol_green = colorRamp2(c(lim_min,0,lim_max), c("#AA4499","white", "#117733"), transparency = 0.01)
diag(est_partcor_plot_subgenes_connected_ordered)=.00001
par(cex=2)
dev.new(); chordDiagramFromMatrix(est_partcor_plot_subgenes_connected_ordered, 
                                  order=c(subgenes.connected.ordered), symmetric=TRUE, 
                                  transparency = 0.1, grid.col="blue", 
                                  col=col_viol_green, annotationTrack=c("grid"), 
                                  keep.diagonal=T, scale=T, reduce=-100, 
                                  preAllocateTracks = list(
                                    track.height = 
                                      max(strwidth(unlist(dimnames(est_partcor_plot_subgenes_connected_ordered)))))) 
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.6)}, bg.border = NA )


