rm(list = ls())
library(caTools)
library(PReMiuM)
library(data.table)
library(RColorBrewer)
library(mcclust)
library(mclust)
library(gtools)

source('simulate_bham_data.R')
source('transform_bham_data.R')

seedySeed <- 12
set.seed(seedySeed)

### Simulating dataset

originalData <- make_fake_data(N_PATIENTS = 6000, N_NULL = 1000)
inputs <- transform_bham_data(originalData$data) # clean and transform before dividing because of empty rows


### split the dataset in two
sample_size <- 0.5 * nrow(inputs$inputData)
indices <- sample(seq_len(nrow(inputs$inputData)), size = sample_size)
train <- inputs$inputData[indices, ]
test <- inputs$inputData[-indices, ]


#Run profile regression on both subsampled datasets
runInfoObjTrain  <- profRegr(yModel="Bernoulli", 
                               xModel="Discrete", 
                               nSweeps=10000, 
                               nClusInit=100, 
                               nBurn=2000,
                               seed = seedySeed,
                               data=train, 
                               output="PremiumOutput/newtrainOutput", 
                               covNames = names(train)[2:(length(names(train)) - 1)], 
                               reportBurnIn = TRUE, 
                               excludeY = TRUE)

runInfoObjTest <- profRegr(yModel="Bernoulli",
                              xModel="Discrete",
                              nSweeps=10000,
                              nClusInit=100,
                              nBurn=2000,
                              seed = seedySeed,
                              data=test,
                              output="PremiumOutput/newtestOutput",
                              covNames = names(test)[2:(length(names(test)) - 1)],
                              reportBurnIn = TRUE,
                              excludeY = TRUE)


zTrain <- fread("PremiumOutput/newtrainOutput_z.txt", header = FALSE)
zMatrixTrain <- as.matrix(zTrain)

zTest <- fread("PremiumOutput/newtestOutput_z.txt", header = FALSE)
zMatrixTest <- as.matrix(zTest)


Rcpp::sourceCpp('makePSM.cpp')

trainPSM <- makePSM(zMatrixTrain)
testPSM <- makePSM(zMatrixTest)

## obtain summary clusterings for both datasets
train$Clusters <- maxpear(trainPSM)$cl
test$Clusters <- maxpear(testPSM)$cl

### How good is the clustering?
print(adjustedRandIndex(train$Clusters, train$group))
print(adjustedRandIndex(test$Clusters, test$group))

#### Maximum likelihood model
clusterParams <- function(data){
  output <- vector(mode = "list")
  mixtureWeights <- vector(mode = "numeric", length = length(unique(data$Clusters)))
  clusterParameters <- matrix(nrow = length(unique(data$Clusters)), ncol = ncol(data) - 3)
  counter <- 1
  for (k in unique(data$Clusters)){
    currentClusterData <- data[data$Clusters == k, 2:(ncol(data)-2)]
    mixtureWeights[counter] <- nrow(currentClusterData)/nrow(data)
    clusterParameters[counter, ] <- colSums(currentClusterData)/nrow(currentClusterData)
    counter <- counter + 1
  }
  output$ClusterParams <- clusterParameters
  output$mixtureWeights <- mixtureWeights
  output
}

clusterParamsTrain <- clusterParams(train)
clusterParamsTest <- clusterParams(test)


clusterPredictions <- function(dataWithLabels, dataWithoutLabels, trainingParams, trainMixWeights){
  posteriorProbs_predictingTest <- matrix(nrow = nrow(dataWithoutLabels), ncol = length(unique(dataWithLabels$Clusters)))
  for(i in 1:nrow(posteriorProbs_predictingTest))
  {
    currentDataForPrediction          <- dataWithoutLabels[i,2:(ncol(dataWithoutLabels)-2)]
    unnormalisedProbs_a                 <- trainingParams[,as.logical(unlist(currentDataForPrediction)), drop = F]
    unnormalisedProbs_b                 <- (1 - trainingParams)[,!as.logical(unlist(currentDataForPrediction)), drop = F]
    unnormalisedProbs                 <- exp(log(trainMixWeights) + rowSums(log(unnormalisedProbs_a)) + rowSums(log(unnormalisedProbs_b)))
    posteriorProbs_predictingTest[i,] <- unnormalisedProbs/sum(unnormalisedProbs)
  }
 apply(posteriorProbs_predictingTest, 1, which.max) 
}



### cluster predictions for both the test and training datasets
predictedClustersForTest <- clusterPredictions(dataWithLabels = train, 
                          dataWithoutLabels = test, 
                          trainingParams = clusterParamsTrain$ClusterParams, 
                          trainMixWeights = clusterParamsTrain$mixtureWeights)

predictedClustersForTrain <- clusterPredictions(dataWithLabels = test, 
                                          dataWithoutLabels = train, 
                                          trainingParams = clusterParamsTest$ClusterParams, 
                                          trainMixWeights = clusterParamsTest$mixtureWeights)


finalTrainClusters <- c(train$Clusters, predictedClustersForTrain)
finalTestClusters <- c(test$Clusters, predictedClustersForTest)

consensusMatrix <- matrix(data = finalTrainClusters, nrow =1)
consensusMatrix <- rbind(consensusMatrix, finalTestClusters)

Rcpp::sourceCpp('makePSM.cpp')
myPSM <- makePSM(consensusMatrix)

heatmapOutput <- pheatmap::pheatmap(myPSM, 
                                    cluster_rows = F, 
                                    cluster_cols = F,
                                    color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), 
                                    cellheight = 0.35, 
                                    cellwidth = 0.35, 
                                    show_rownames = F)

