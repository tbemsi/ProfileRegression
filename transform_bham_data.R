semi.deterministic <- function(probs=c(0.7,0.3), size=10){
  x1 <- sample(c(0,1), size/2, replace=TRUE, prob=probs)
  x2 <- sample(c(0,1), size/2, replace=TRUE, prob=rev(probs))
  c(x1, x2)
}

continuous <- function(size=10, means = 1, distance = 0.5){
  x1 <- rnorm(size/4, mean = means + 0*distance)
  x2 <- rnorm(size/4, mean = means + 1*distance)
  x3 <- rnorm(size/4, mean = means + 2*distance)
  x4 <- rnorm(size/4, mean = means + 3*distance)
  c(x1, x2, x3, x4)
}

transform_bham_data <- function(inputs){
  ready_data <- vector(mode = "list")
  ready_data$inputData <- inputs[1:27]
  random.outcome <- sample(c(0,1), nrow(inputs), replace=TRUE, prob=c(0.5,0.5))
  deterministic.outcome <- c(rep(0, nrow(inputs)/2), rep(1, nrow(inputs)/2))
  semi.deterministic.outcome.one <- semi.deterministic(probs = c(0.6, 0.4), size=nrow(inputs))
  semi.deterministic.outcome.two <- semi.deterministic(probs = c(0.7, 0.3), size=nrow(inputs))
  semi.deterministic.outcome.three <- semi.deterministic(probs = c(0.8, 0.2), size=nrow(inputs))
  semi.deterministic.outcome.four <- semi.deterministic(probs = c(0.9, 0.1), size=nrow(inputs))
  continuous.outcome.one <- continuous(nrow, means=1, distance = 1)
  continuous.outcome.two <- continuous(nrow, means=1, distance = 3)
  continuous.outcome.three <- continuous(nrow, means=1, distance = 5)
  continuous.outcome.four <- continuous(nrow, means=1, distance = 1)
  covnames <- c("group", rep(0, 26))
  for (i in 2:27){
    covnames[i] <- paste("Variable",i, sep="")
  }
  ready_data$covNames <- covnames
  colnames(ready_data$inputData) <- covnames
  ready_data$inputData$RandomOutcome <- random.outcome
  ready_data$inputData$DetOutcomeOne <- deterministic.outcome
  ready_data$inputData$DetOutcomeTwo <- semi.deterministic.outcome.one
  ready_data$inputData$DetOutcomeThree <- semi.deterministic.outcome.two
  ready_data$inputData$DetOutcomeFour <- semi.deterministic.outcome.three
  ready_data$inputData$DetOutcomeFive <- semi.deterministic.outcome.four
  ready_data$inputData$ContinuousOne <- continuous.outcome.one
  ready_data$inputData$ContinuousTwo <- continuous.outcome.two
  ready_data$inputData$ContinuousThree <- continuous.outcome.three
  ready_data$inputData$ContinuousFour <- continuous.outcome.four
  ready_data$xModel <- "Discrete"
  ready_data$yModel <- "Bernoulli"
  ready_data$nCovariates <- 26 
  ready_data
}