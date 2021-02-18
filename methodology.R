rm(list = ls())
library(caTools)
library(PReMiuM)
library(data.table)
library(RColorBrewer)
library(mcclust)
library(mclust)
library(gtools)

source('~/Documents/BSU/code/simulate_bham_data.R')
source('~/Documents/BSU/code/transform_bham_data.R')
setwd('~/Documents/BSU/code/')
seedy.seed <- 12
set.seed(seedy.seed)

### Simulating dataset

df <- make_fake_data(N_PATIENTS = 6000, N_NULL = 1000)
inputs <- transform_bham_data(df$data) # clean and transform before dividing because of missing rows


### split the dataset in two
sample_size <- 0.5 * nrow(inputs$inputData)
indices <- sample(seq_len(nrow(inputs$inputData)), size = sample_size)
train <- inputs$inputData[indices, ]
test <- inputs$inputData[-indices, ]


### transform each subsampled dataset
data_transformation <- function(inputs){
  ready_data <- vector(mode = "list")
  rownames(inputs) <- seq(nrow(inputs))
  ready_data$inputData <- inputs
  covnames <- c("group", rep(0, 26))
  for (i in 2:27){
    covnames[i] <- paste("Variable",i, sep="")
  }
  ready_data$covNames <- covnames
  colnames(ready_data$inputData) <- covnames
  ready_data$inputData$Outcome <- inputs$Outcome
  ready_data$xModel <- "Discrete"
  ready_data$yModel <- "Bernoulli"
  ready_data$nCovariates <- 26 
  ready_data
}


new.train <- data_transformation(train)
new.test <- data_transformation(test)

#Run profile regression on both subsampled datasets
runInfoObj.train   <- profRegr(yModel=new.train$yModel, 
                                    xModel=new.train$xModel, 
                                    nSweeps=10000, 
                                    nClusInit=100, 
                                    nBurn=2000,
                                    seed = seedy.seed,
                                    data=new.train$inputData, 
                                    output="newtrainOutput", 
                                    covNames = new.train$covNames, 
                                    reportBurnIn = TRUE, 
                                    excludeY = TRUE)

z.train <- fread("newtrainOutput_z.txt", header = FALSE)
zMatrix <- as.matrix(z.train)


Rcpp::sourceCpp('makePSM.cpp')
myPSM.train <- makePSM(zMatrix)


set.seed(seedy.seed)
myClustersReduced.train <- maxpear(myPSM.train)$cl
train$Clusters <- myClustersReduced.train

#### Maximum likelihood model

compute_theta <- function(df){
  number.clusters <- length(unique(df$Clusters))
  theta <- matrix(nrow = number.clusters, ncol = 26)
  for (i in 1:number.clusters){
    cluster.weight <- length(which(df$Clusters ==i))
    for (j in 2:27){
      theta[i, j-1] <- 1/cluster.weight * sum(subset(train, Clusters==i)[, j])
    }
  }
  theta
}


compute_likelihood <- function(observation, cluster, theta){
  likelihood = 1
  for (j in 1:length(observation)){
    likelihood = likelihood * theta[cluster, j]**(observation[j]) * (1 - theta[cluster,j])**(1 - observation[j])
  }
  likelihood
}

posterior_prob <- function(observation, df, theta){
  number.clusters <- length(unique(df$Clusters))
  cluster.weights <- as.data.frame(table(df$Clusters))
  posterior.probs <- compute_likelihood(observation, 1:number.clusters, theta)*as.numeric(cluster.weights$Freq)
  which.max(posterior.probs)
}



###### Testing data set

runInfoObj.test   <- profRegr(yModel=new.test$yModel,
                         xModel=new.test$xModel,
                         nSweeps=10000,
                         nClusInit=100,
                         nBurn=2000,
                         data=new.test$inputData,
                         output="newtestOutput",
                         covNames = new.test$covNames,
                         reportBurnIn = TRUE,
                         excludeY = TRUE)

z.test <- fread("newtestOutput_z.txt", header = FALSE)
zMatrix <- as.matrix(z.test)


Rcpp::sourceCpp('makePSM.cpp')
myPSM.test <- makePSM(zMatrix)

set.seed(seedy.seed)
myClustersReduced <- maxpear(myPSM.test)$cl
number.clusters <- length(unique(myClustersReduced))

test$Clusters <- myClustersReduced


########----------------------------------------------

theta.train <- compute_theta(train)
theta.test <- compute_theta(test)
observation.train = as.numeric(train[1, ][2:27])
observation.test = as.numeric(test[1, ][2:27])
posterior_prob(observation.train, train, theta.train)
posterior_prob(observation.test, test, theta.test)


#### take what has been trained on the training set, and apply to the testing set
sth.test <- rep(0, nrow(test))
for (i in 1:nrow(train)){
  observation <- as.numeric(test[i, ][2:27])
  sth.test[i] <- posterior_prob(observation, train, theta.train)
}
new.test <- test[, c(1:28)]
new.test$Clusters <- sth.test

final.df.one <- rbind(train, new.test)


sth.train <- rep(0, nrow(train))
for (i in 1:nrow(test)){
  observation <- as.numeric(train[i, ][2:27])
  sth.train[i] <- posterior_prob(observation, test, theta.test)
}
new.train <- train[, c(1:28)]
new.train$Clusters <- sth.train

final.df.two <- rbind(test, new.train)


ARI <- adjustedRandIndex(final.df.one, final.df.two)


##### Produce matrix. Fit test using the training data, get the clustering allocations
##### Matrix[i,j] = sum(Indicator(person is allocated in cluster i from premium, and cluster j from MLE))
##### Maybe also use ARIs between different clusterings

#### New likelihood function which takes into account the outcome















