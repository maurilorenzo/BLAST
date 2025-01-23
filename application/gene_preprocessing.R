# R SCRIPT TO POST-PROCESS DATA USED IN APPLICATION IN SECTION 5 OF THE PAPER

genedata=readRDS('application/data/genedata.rds')
# post processing from https://github.com/noirritchandra/SUFA/blob/main/vignettes/Genedata_application.R
set.seed(123)

bulk.data= genedata$bulk
array.data= genedata$array1
array.data2= genedata$array2

Ident.bulk= genedata$bulk.types
Ident.array= genedata$array1.types
Ident.array2= genedata$array2.types

dat=ExpressionSet(assayData = bulk.data)
filter.dat=varFilter(dat, var.cutoff = cutoff)
hvgs_in_bulk=rownames(filter.dat)
rm(filter.dat,dat)

dat=ExpressionSet(assayData = array.data)
filter.dat=varFilter(dat, var.cutoff = cutoff)
hvgs_in_array=rownames(filter.dat)
rm(filter.dat,dat)

dat=ExpressionSet(assayData = array.data2)
filter.dat=varFilter(dat, var.cutoff = cutoff)
hvgs_in_array2=rownames(filter.dat)
rm(filter.dat,dat)

genes.use=Reduce(intersect,list(hvgs_in_array,hvgs_in_array2,hvgs_in_bulk))
length(genes.use)

Y=list(t(array.data[genes.use,]), t(array.data2[genes.use,]),t(bulk.data[genes.use,] ))
Ident.list=list(Ident.array,Ident.array2,Ident.bulk)

Y_cent=mapply(function(Idents,y){
  cell.types=unique(Idents)
  y_cent=lapply(cell.types, function(ident,Identt,Y){
    scale(Y[which(Identt==ident),],scale=F)
  } ,Y=y,Identt=Idents)
  Reduce(rbind,y_cent)
}, Ident.list, Y )

sd.median=median( sapply(Y_cent, function(y) apply(y,2,sd)))
Y_cent_scaled=lapply(Y_cent, "/",   sd.median)

Y_1 <- Y_cent_scaled[[1]]; dim(Y_1)
Y_2 <- Y_cent_scaled[[2]]; dim(Y_2)
Y_3 <- Y_cent_scaled[[3]]; dim(Y_3)

Y <- list(Y_1, Y_2, Y_3)

n_1 <- nrow(Y_1)
n_2 <- nrow(Y_2)
n_3 <- nrow(Y_3)

set.seed(123)
test_gene_1 <- sample(1:n_1, floor(0.2*n_1))
test_gene_2 <- sample(1:n_2, floor(0.2*n_2))
test_gene_3 <- sample(1:n_3, floor(0.2*n_3))

Y_1_train_gene <- as.matrix(Y_1[-test_gene_1,])
Y_2_train_gene <- as.matrix(Y_2[-test_gene_2,])
Y_3_train_gene <- as.matrix(Y_3[-test_gene_3,])
Y_train_gene <- list(Y_1_train_gene, Y_2_train_gene, Y_3_train_gene )  

Y_1_test_gene <- as.matrix(Y_1[test_gene_1,])
Y_2_test_gene <- as.matrix(Y_2[test_gene_2,])
Y_3_test_gene <- as.matrix(Y_3[test_gene_3,])

Y_test_gene <- list(Y_1_test_gene, Y_2_test_gene, Y_3_test_gene)  

p <- ncol(Y_1_test_gene)
set.seed(123)
impute_factor_index_gene <- sample(1:p, floor(0.5*p))
p_test <- p-length(impute_factor_index_gene)

