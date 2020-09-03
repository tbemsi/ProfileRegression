library('gtools')
library('Rlab')

simulate_data <- function(number_conditions=3,
                          number_clusters = 3,
                          number_patients = 5,
                          alpha = 0.1
){
  conditions <- letters[1:number_conditions] ### we will just call the conditions a, b, c,...
  cluster_names <- LETTERS[1:number_clusters]
  clusters <- sample(cluster_names, number_patients, replace=T) ### This tells me what cluster each person belongs to
  
  ### sample according to different probabilities of belonging to different clusters
  ### look at premium code for simulating data
  
  ### work out how to modify parameters of clusSummaryBernouilliDiscrete() function
  
  ### write new version of clus... fn which gives binary data (nCategories = 2), no fixed effects (nFixedEffects = 0)
  
  probability_matrix <- rdirichlet(number_clusters, rep(alpha, number_conditions)) ### generate vector of probabilities for each cluster
  rownames(probability_matrix) <- cluster_names
  simulated_data <- matrix(0, number_patients, number_conditions)
  
  for (i in 1:number_patients){
    random_vector <- probability_matrix[clusters[i],]
    simulated_data[i,] <- rbern(number_conditions, prob=random_vector)
    colnames(simulated_data) <- conditions
  }
  simulated_data <- cbind(clusters, as.data.frame(simulated_data))
  simulated_data
}

simulate_data(number_conditions = 10, number_clusters = 3, number_patients = 30)

clusSummaryBernoulliDiscrete_mine <- function(outcomeType = "Bernouilli",
                                              covariateType = "Discrete",
                                              nCovariates = 5,
                                              nCategories = rep(2, nCovariates),
                                              nFixedEffects = 0,)