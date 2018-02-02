consensusNet <- function(data, organism="hsapiens",bootstrapNum=100, naPer=0.5, meanPer=0.8,varPer=0.8,method="rank_unsig",value=3/1000,pth=1e-6, 
                         nMatNet=2, nThreads=4){
# there are two levels of parallelism here
# level 1: the number of concurrent running MatNet processes
# level 2: for each MatNet process, how many threads can be use 


#  entrezId <- IDMapping(organism=organism,dataType="list",
#                       inputGene=rownames(data), sourceIdType="genesymbol", 
#                       targetIdType="entrezgene", hostName=hostName)
#
#  entrezId <- entrezId$mapped
#  data <- data[entrezId[,1],]
#  rownames(data) <- entrezId[,4]

  cl <- makeCluster(nMatNet, outfile="")
  registerDoParallel(cl)

  sampleName <- colnames(data)
  
  message("Creating random networks.....")
  bN <- foreach(i=1:bootstrapNum) %dopar% {
    ransample <- sample(sampleName,length(sampleName),replace=TRUE)
    ranData <- data[,ransample]

    if(method=="rank"){
      ranNetwork <- MatNet(inputMat=ranData, naPer=naPer, meanPer=meanPer,varPer=varPer,corrType="spearman",matNetMethod=method,rankBest=value,nThreads=nThreads)
    }

    if(method=="value"){
      ranNetwork <- MatNet(inputMat=ranData, naPer=naPer, meanPer=meanPer,varPer=varPer,corrType="spearman",matNetMethod=method,valueThr=value,nThreads=nThreads)
    }

    if(method=="rank_unsig"){
      ranNetwork <- MatNet(inputMat=ranData, naPer=naPer, meanPer=meanPer,varPer=varPer,corrType="spearman",matNetMethod="rank",rankBest=value,netFDRThr=1,nThreads=nThreads)
    }
    return(ranNetwork)
  }

  stopCluster(cl)

  message("Calculating consensus network .....")
  consensusNetwork <- .calculateConsensusNetwork(ranNetList=bN,pth=pth)
  return(consensusNetwork)
}


.normalizeNetwork <- function(network){
    x <- network[,c(2,1)]
    colnames(x) <- colnames(network)
    x <- rbind(x,network)
    x <- x[x[,1]<x[,2],]
    x <- unique(x)
    return(x)
}


.calculateConsensusNetwork <- function(ranNetList,pth){
    bootStrapNum <- length(ranNetList)
    ranEdgeNum <- c()
    
    allEdge <- c()
    
    for(i in c(1:bootStrapNum)){
        ranNet <- ranNetList[[i]]
        ranNet <- .normalizeNetwork(ranNet)
        ranEdgeNum <- c(ranEdgeNum,nrow(ranNet))
        
        ranNet <- paste(ranNet[,1],ranNet[,2],sep="_")
        allEdge <- c(allEdge,ranNet)
    }
    
    totalEdgeNum <- length(unique(allEdge))
    
    cat("Bonferroni 0.05:",0.05/totalEdgeNum,"\n")
    
    prob <- ranEdgeNum/totalEdgeNum
    
    mu <- sum(prob)
    sigma <- sqrt(sum(prob*(1-prob)))
    
    edgeRe <- table(allEdge)
    edgeReZ <- (edgeRe-mu)/sigma
    edgeReP <- pnorm(edgeReZ,mean=0,sd=1,lower.tail=FALSE)
    
    edgeReP <- edgeReP[edgeReP<pth]
    if(length(edgeReP)==0){
        return(NULL)
    }else{
        consensusNetwork <- names(edgeReP)
        consensusNetwork <- do.call(rbind,strsplit(consensusNetwork,"_"))
        return(consensusNetwork)
    }
}
