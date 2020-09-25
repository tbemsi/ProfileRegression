generate_cluster_data <- function(nClusters = 3,
                                  nCovariates = 15,
                                  probs_on = c(0.9,0.1),
                                  probs_off = c(0.1,0.9),
                                  theta = rep(0, nClusters)){
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
    clusterData[[i]] <- list(theta=theta[i], covariateProbs = covariateProbs)
    
  }
  clusterData
}

clusSummaryBernoulliDiscrete_mine <- function(outcomeType = "Bernoulli",
                                              covariateType = "Discrete",
                                              nCovariates = 5,
                                              nClusters = 3,
                                              nCategories = rep(2, nCovariates),
                                              nFixedEffects = 0,
                                              clusterSizes = c(200, 200, 200),
                                              probs_on = c(0.8,0.2),
                                              probs_off = c(0.2,0.8)){
  list(outcomeType = outcomeType, covariateType = covariateType, 
       nCovariates = nCovariates, nCategories = nCategories, nFixedEffects = nFixedEffects, 
       fixedEffectsCoeffs = c(0.1, -0.5), missingDataProb = 0, 
       nClusters = nClusters, clusterSizes = clusterSizes, 
       includeCAR = FALSE, TauCAR = 100, clusterData = generate_cluster_data(nClusters = nClusters,
                                                                             nCovariates = nCovariates,
                                                                             probs_on = probs_on,
                                                                             probs_off = probs_off)) 
}
