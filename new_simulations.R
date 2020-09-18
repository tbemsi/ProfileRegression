rm(list = ls())
graphics.off()
setwd('/Users/bmsbm/Documents/BSU/MM') 
library(PReMiuM)
library(data.table)
library(RColorBrewer)
library(mcclust)
library(gtools)
set.seed(2)

source('clusSummaryBernouilli_new.R')
source('generate_psm.R')

#-----------------------1. Setting cluster sizes manually -----------------------------------------------
### generate some data
nCovariates <- 5
nClusters <- 3
ClusterSizes <- c(300, 100, 100)
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete_mine(nCovariates = nCovariates,
                                                                   nClusters = nClusters,
                                                                   clusterSizes = ClusterSizes))
cluster_labels <- generate_cluster_labels(ClusterSizes)
annot <- data.frame(Group = as.factor(cluster_labels)) 
rownames(annot) <- rownames(inputs$inputData) <- seq(1, nrow(inputs$inputData))

### Visualise the data:
pheatmap::pheatmap(inputs$inputData[,2:ncol(inputs$inputData)],  
                   cluster_rows = F,
                   cluster_cols = F,
                   col = c("white", "black"),
                   annotation_row = annot)


### Apply Premium on generated data
#generate_psm(inputs, filename = "psm_manual_weights")

myOrder2 <- generate_psm(inputs, filename = "psm_manual_weights") # this is the order of the PSM


### now I want to do a heatmap of the data with myOrder2
sth <- inputs$inputData[myOrder2[1], ]
for (i in 2:length(myOrder2)){
  sth <- rbind(sth, inputs$inputData[myOrder2[i], ])
}
ordered_df <- as.data.frame(sth)

myPSMHeatmap <- pheatmap::pheatmap(ordered_df,  
                                   cluster_rows = F,
                                   cluster_cols = F,
                                   col = c("white", "black"),
                                   annotation_row = annot,
                                   filename = "ordered_df_manual_weights.pdf")


#2. -------------------- Using dirichlet process to generate cluster weights -------------------
nClusters <- 3
nPatients <- 600
alphaParameter <- 1
mixtureWeights <- rdirichlet(1, rep(alphaParameter,nClusters))
clusterSizes <- rmultinom(1, nPatients, mixtureWeights)  
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete_mine(nCovariates = nCovariates,
                                                                   nClusters = nClusters,
                                                                   clusterSizes = ClusterSizes))

cluster_labels <- generate_cluster_labels(ClusterSizes)
annot <- data.frame(Group = as.factor(cluster_labels)) 
rownames(annot) <- rownames(inputs$inputData) <- seq(1, nrow(inputs$inputData))

### Visualise the data:
pheatmap::pheatmap(inputs$inputData[,2:ncol(inputs$inputData)],  
                   cluster_rows = F,
                   cluster_cols = F,
                   col = c("white", "black"),
                   annotation_row = annot)


### Apply Premium on generated data
#generate_psm(inputs)

myOrder2 <- generate_psm(inputs, filename="psm_dirichlet") # this is the order of the PSM


### now I want to do a heatmap of the data with myOrder2
sth <- inputs$inputData[myOrder2[1], ]
for (i in 2:length(myOrder2)){
  sth <- rbind(sth, inputs$inputData[myOrder2[i], ])
}
ordered_df <- as.data.frame(sth)

myPSMHeatmap <- pheatmap::pheatmap(ordered_df,  
                                   cluster_rows = F,
                                   cluster_cols = F,
                                   col = c("white", "black"),
                                   annotation_row = annot,
                                   filename = "ordered_df_dirichlet.pdf")


#3.----------------------- Traceplot of alpha parameter ----------------------------------------
alphas <- read.delim("myOutput_unsupervised_alpha.txt", header = FALSE, sep = "\t")
alphas <- alphas[seq(1, nrow(alphas), by=50),] # thinning the alphas
plot(alphas, type ='l', ylab='alpha')


# increase dimensionality of data and see 10 covariates, maybe 1000 people
# calculated Adjusted Rand Index (mcclust package) between final clustering and true clustering
# distribution of ARIs - each row of z matrix and the ground truth - maybe boxplots?









