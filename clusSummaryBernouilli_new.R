rm(list = ls())
graphics.off()
#setwd('/Users/bmsbm/Documents/BSU/MM')
library(PReMiuM)
library(data.table)
library(RColorBrewer)
library(mcclust)
set.seed(2)


generate_cluster_data <- function(nClusters = 3,
                                  nCovariates = 5,
                                  probs_on = c(0.9,0.1),
                                  probs_off = c(0.1,0.9)){
  clusterData <- list()
  for (i in 1:nClusters){
    number_switched_on <- 3 #sample(1:nCovariates, 1) ### deciding how many variables to be switched on in cluster i
    covariateProbs <- vector(mode = "list", length = nCovariates)
    on_positions <- sample(1:nCovariates, number_switched_on) # positions of 'on' indices in covariateProbs
    off_positions <- setdiff(1:nCovariates, on_positions) # positions of 'off' indices in covariateProbs
    for (j in on_positions){
      covariateProbs[[j]] <- probs_on
    }
    for (k in off_positions){
      covariateProbs[[k]] <- probs_off
    }
    clusterData[[i]] <- list(theta=0, covariateProbs = covariateProbs)
    
  }
  clusterData
}

clusSummaryBernoulliDiscrete_mine <- function(outcomeType = "Bernouilli",
                                              covariateType = "Discrete",
                                              nCovariates = 5,
                                              nClusters = 3,
                                              nCategories = rep(2, nCovariates),
                                              nFixedEffects = 0,
                                              probs_on = c(0.8,0.2),
                                              probs_off = c(0.2,0.8)){
  list(outcomeType = outcomeType, covariateType = covariateType, 
       nCovariates = nCovariates, nCategories = nCategories, nFixedEffects = nFixedEffects, 
       fixedEffectsCoeffs = c(0.1, -0.5), missingDataProb = 0, 
       nClusters = nClusters, clusterSizes = c(200, 200, 200), 
       includeCAR = FALSE, TauCAR = 100, clusterData = generate_cluster_data(nClusters = nClusters,
                                                                             nCovariates = nCovariates,
                                                                             probs_on = probs_on,
                                                                             probs_off = probs_off)) 
}

### generate some data
inputs <- generateSampleDataFile(clusSummaryBernoulliDiscrete_mine())

myData <- inputs$inputData
myDataInput <- myData[,tail(colnames(myData), 5)] # Remember to change this to nCovariates when the fxn is working

#Visualise the data:
pheatmap::pheatmap(myDataInput, 
                   cluster_rows = F, 
                   cluster_cols = F,
                   col = c("white", "black"),
                   #annotation_row = annot, 
                   cellheight = 0.35, 
                   show_rownames = F)


### Apply Premium on generated data


### **Goal**: generate some data, visualize it including the true cluster labels as an annotation.
### Applying Premium on the generated data to see how it performs

### Merge my notes and Paul's

