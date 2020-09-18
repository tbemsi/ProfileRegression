generate_cluster_labels <- function(ClusterSizes){
  cluster_labels <- rep(1, ClusterSizes[1])
  for (i in 2:length(ClusterSizes)){
    cluster_labels <- c(cluster_labels,(rep(i,ClusterSizes[i])))
  }
  cluster_labels
}

generate_psm <- function(inputs, filename = 'psm.pdf'){
  runInfoObj <- profRegr(yModel = inputs$yModel,
                         xModel = inputs$xModel,
                         nSweeps = 100000,
                         nClusInit=100,
                         nBurn=2000,
                         data=inputs$inputData[,2:ncol(inputs$inputData)],
                         output="myOutput_unsupervised",
                         covNames = names(inputs$inputData[,2:ncol(inputs$inputData)]),
                         reportBurnIn = TRUE, 
                         excludeY = TRUE)
  
  # trace plot to decide how many iterations count as burn-in
  globalParsTrace(runInfoObj,
                  parameters = "nClusters",
                  plotBurnIn=FALSE,
                  whichBeta=1) 
  
  z <- fread("myOutput_unsupervised_z.txt", header = FALSE)
  zMatrix <- as.matrix(z)
  zMatrix <- zMatrix[ceiling(nrow(zMatrix)/2):nrow(zMatrix),]
  
  
  #Thin the MCMC output by taking only every 500-th draw:
  zMatrix <- zMatrix[seq(1, nrow(zMatrix), by=500),]
  
  
  ### Calculate the PSM:
  Rcpp::sourceCpp('makePSM.cpp')
  myPSM <- makePSM(zMatrix)
  save(myPSM, file = "myOutput_PSM.RDa")
  
  set.seed(1)
  myData <- inputs$inputData
  myDataReduced     <- inputs$inputData[,2:ncol(inputs$inputData)]
  myClustersReduced <- maxpear(myPSM)$cl
  
  
  # Create annotation for heatmaps
  annot <- data.frame(Cluster = as.factor(myClustersReduced))
  rownames(annot) <- rownames(myPSM) <- rownames(myDataReduced)
  plottingOrder <- order(myClustersReduced)
  graphics.off()
  
  
  # Plot the PSM (with rows and columns ordered according to cluster label), together with the cluster labels  
  pdf(file = filename)
  heatmapOutput <- pheatmap::pheatmap(myPSM[plottingOrder,plottingOrder], 
                                      cluster_rows = T, 
                                      cluster_cols = T,
                                      color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100), 
                                      annotation_row = annot, 
                                      cellheight = 0.35, 
                                      cellwidth = 0.35, 
                                      show_rownames = F)
  dev.off()
  myOrder <- heatmapOutput$tree_row$order
  myOrder
}
