make_fake_data <- function(
  K=3, 
  allow_overlap=TRUE, 
  balanced=TRUE, 
  N_PATIENTS=10000, 
  N_NULL=1000,
  NOISE_PARAMETER = -Inf, # No noise
  MEAN_VAL = -0.7, # within cluster mean value on qnorm scale
  COR_VAL = 0.7 # within cluster correlation value
){
  
  if(balanced==FALSE){
    # imbalanced scenario
    balance_type = "unbalanced"
    patient_groups <- sample(LETTERS[1:K], N_PATIENTS, prob=rexp(K), replace=T)  #1.
  } else {
    # balanced scenario - uncomment as required
    balance_type = "balanced"
    patient_groups <- sample(LETTERS[1:K], N_PATIENTS, replace=T)
  }
  
  # for each of these groups, allocate a cluster
  # assign a random sample of letters to each group member, with size determined by pmax(2, rpois(n=1, lambda=5))
  
  disease_groups <- list()
  
  if(allow_overlap){
    # allow overlap in this case
    for(i in 1:K) disease_groups[[i]] <- sort(sample(letters, pmax(2, rpois(n=1, lambda=5))))
  } else {
    # do not allow overlap in this case
    to_sample <- letters
    for(i in 1:K) {
      disease_groups[[i]] <- sort(sample(to_sample, pmax(2, rpois(n=1, lambda=5))))
      to_sample <- setdiff(to_sample, disease_groups[[i]])
    }
  }
  
  names(disease_groups) <- LETTERS[1:K]
  
  # now we've got a list of which disease clusters go with each patient group, simulate this...
  # 25 different disease groups comprising 26 different diseases.
  # let's use the probit model to generate presence/absence.
  # assign a mean value for the disease if present and then a correlation matrix.
  
  require(mvtnorm)
  fake_output <- NULL
  # add a background probability for each condition
  for(disease_group_name in names(disease_groups)){
    disease_group <- disease_groups[[disease_group_name]]
    # how many patients have the disease within the group
    N <- sum(patient_groups == disease_group_name)
    n <- length(disease_group)
    # sample a within-group correlation value - we will say this is constant within the group
    # cor_val <- runif(1, 0.3, 1) # 0.7
    
    # Tom modification for deterministic CORRELATION VALUE
    #cor_val <- 0.7
    cor_val <- COR_VAL
    
    mean_vals <- rep(-Inf, 26)  ### 3
    w <- which(letters %in% disease_group)
    
    # mean_vals[w] <- pmin(rnorm(n, -0.7, 0.7), -0.5) # value of the probit - set this with a cap at 30%
    
    # Tom modification for deterministic MEAN PREVALENCE
    #mean_vals[w] <- -0.7
    mean_vals[w] <- MEAN_VAL
    
    # form the correlation matrix
    covmat <- diag(26)
    covmat[w, w] <- matrix(cor_val, n, n)
    diag(covmat) <- 1
    
    tmp <- rmvnorm(N, mean_vals, covmat)  ## 2
    tmp <- 1*(pnorm(tmp) > runif(length(tmp)))
    tmp <- data.frame(disease_group_name, tmp)
    colnames(tmp) <- c("group", letters)
    fake_output <- rbind(fake_output, tmp)
  }
  
  # now add say N_NULL zeroes for group Z
  null_df <- data.frame("Z", matrix(0, N_NULL, 26))
  colnames(null_df) <- c("group", letters)
  fake_output <- rbind(fake_output, null_df)
  
  # now add background noise - set the background noise level to a proportion of the mean value and add uncorrelated noise
  # adj - offset, if this is lower the noise will be lower
  make_noise <- function(adj){
    #print(adj)
    mean_vals <- colMeans(fake_output[,-1])
    mean_vals[] <- adj
    covmat <- diag(26)
    tmp <- rmvnorm(nrow(fake_output), mean_vals, covmat)
    tmp <- 1*(pnorm(tmp) > runif(length(tmp)))
    #print(head(tmp))
    tmp
  }
  
  
  ## low noise
  # retain group identifier so that conditions in each cluster can be added
  fake_plus_noise <- data.frame(group=fake_output[ , 1], pmax(make_noise(NOISE_PARAMETER), as.matrix(fake_output[ , -1])))
  # add conditions in cluster to each patient record
  fake_plus_noise <- transform(fake_plus_noise, conditions=sapply(disease_groups, paste0, collapse=' ')[group])
  fake_plus_noise$conditions <- as.character(fake_plus_noise$conditions)
  
  list(data=fake_plus_noise, disease_groups=disease_groups)
  #list(data=fake_plus_noise, disease_groups=disease_groups)
}


#zz1 <- make_fake_data(COR_VAL=1.0)

#image(as.matrix(zz1$data[,letters]))

#zz2 <- zz1$data
#table(zz2$conditions)
