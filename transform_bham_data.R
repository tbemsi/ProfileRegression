transform_bham_data <- function(inputs, outcomeProbs = c(0.5, 0.5)){
  ready_data <- vector(mode = "list")
  inputs$totals <- apply(inputs[, c(2:27)], 1, sum)
  inputs <- subset(inputs, totals > 0)
  rownames(inputs) <- seq(nrow(inputs))
  ready_data$inputData <- inputs[1:27]
  outcome <- c(sample(c(0,1), floor(nrow(inputs)/2), replace=TRUE, prob=outcomeProbs),
               sample(c(0,1), ceiling(nrow(inputs)/2), replace=TRUE, prob=rev(outcomeProbs)))
  covnames <- c("group", rep(0, 26))
  for (i in 2:27){
    covnames[i] <- paste("Variable",i-1, sep="")
  }
  ready_data$covNames <- covnames
  colnames(ready_data$inputData) <- covnames
  ready_data$inputData$Outcome <- outcome
  ready_data$xModel <- "Discrete"
  ready_data$yModel <- "Bernoulli"
  ready_data$nCovariates <- 26 
  ready_data
}

