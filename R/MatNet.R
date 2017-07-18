MatNet <-
function(inputMat, collapse_mode="maxSD", naPer=0.7, meanPer=0.8, varPer=0.8, corrType="spearman", matNetMethod="rank", valueThr=0.5, rankBest=0.003, networkType="signed", netFDRMethod="BH", netFDRThr=0.05, idNumThr=(-1), nThreads=3){
	
    #require(igraph) || stop("Package igraph version 0.6 is required!")
    #require(WGCNA) || stop("Package WGCNA version 1.34.0 is required!")
    #require(tools) || stop("Package tools version 3.0.0 is required!")
    
    #load expression data
    
    re <- testFileFormat(inputMat=inputMat,collapse_mode=collapse_mode)
    inputMat <- re$inputMat
    inputMat <- as.matrix(inputMat)
    
    if(naPer>1 || naPer<0 || class(meanPer)!="numeric"){
        stop("The input 'naPer' is invalid! 'naPer' should be a positive number and less than 1!\n")
    }
    
    if(meanPer>1 || meanPer<0 || class(meanPer)!="numeric"){
        stop("The input 'meanPer' is invalid! 'meanPer' should be a positive number and less than 1!\n")
    }
    
    if(varPer>1 || varPer<0 || class(varPer)!="numeric"){
        stop("The input 'varPer' is invalid! 'varPer' should be a positive number and less than 1!\n")
    }
    
    
    if(length(which(corrType %in% c("pearson","spearman")))==0){
		stop("The input 'corrType' is invalid! Please select an option from 'pearson' and 'spearman'!\n")
	}
    
    if(length(which(matNetMethod %in% c("value","rank","directed")))==0){
		stop("The input 'matNetMethod' is invalid! Please select an option from 'value','rank' and 'directed'!\n")
	}
    
    if(valueThr>1 || valueThr<0 || class(valueThr)!="numeric"){
        stop("The input 'valueThr' is invalid! 'valueThr' should be a positive number and less than 1!\n")
    }
    
    if(rankBest<0 || rankBest>1){
        stop("The input 'rankBest' is invalid! 'rankBest' should be from 0 to 1!\n")
    }
    
    if(length(which(networkType %in% c("signed","unsigned")))==0){
		stop("The input 'networkType' is invalid! Please select from 'signed' and 'unsigned'!\n")
	}
    
    if(length(which(netFDRMethod %in% c("holm","hochberg","hommel","bonferroni","BH","BY","none")))==0){
		stop("The input 'netFDRMethod' is invalid! Please select a method from 'holm','hochberg', 'hommel', 'bonferroni', 'BH', 'BY' and 'none'!\n")
	}
    
    if(netFDRThr>1 || class(netFDRThr)!="numeric"){
        stop("The input 'netFDRThr' is invalid! 'netFDRThr' should be a number and less than 1!\n")
    }
    
    if(idNumThr<(-1) || class(idNumThr)!="numeric"){
        stop("The input 'idNumThr' is invalid! 'idNumThr' should be -1 or a positive number!\n")
    }
    
    if(nThreads<1 || class(nThreads)!="numeric"){
        stop("The input 'nThreads' is invalid! 'nThreads' should be a number and greater than 1!\n")
    }
    
    missing <- 0
    
	if(length(which(is.na(inputMat)))>0){
		cat("The input data contain ",nrow(inputMat)," ids and ",ncol(inputMat)," samples. The data contain missing values.\n",sep="")
		missing <- 1
	}else{
		cat("The input data contain ",nrow(inputMat)," ids and ",ncol(inputMat)," samples. No missing values in the data.\n",sep="")
	}
    
    #filter genes
    #cat("Remove low expressed genes and less variable genes...\n")
    
    if(missing==1){
    
        naG <- apply(inputMat,1,.calculateNANum)
        if(length(naG[naG>naPer])>0){
            cat(paste("Each of ",length(naG[naG>naPer])," ids has over ",naPer*100,"% missing values in all ",ncol(inputMat)," samples.The function will remove these ids.\n",sep=""))
            naG <- naG[naG<=naPer]
            filterData <- inputMat[names(naG),]
        }else{
            filterData <- inputMat
        }
    }else{
        filterData <- inputMat
    }
    
    meanG <- apply(filterData,1,mean,na.rm=T)
    sdG <- apply(filterData,1,sd,na.rm=T)
    
    if(length(which(sdG==0))>0){
        x <- sdG[sdG==0]
        cat(paste("The function removed ",length(x)," ids with standard deviation 0.\n",sep=""))
        sdG <- sdG[sdG!=0]
        meanG <- meanG[names(sdG)]
        filterData <- filterData[names(sdG),]
    }
    
    
    meanG <- sort(meanG,decreasing=T)
    meanG <- meanG[c(1:round(meanPer*nrow(filterData)))]
    sdG <- sdG[names(meanG)]
    sdG <- sort(sdG,decreasing=T)
    sdG <- sdG[c(1:round(varPer*length(meanG)))]
    
    if(idNumThr!=(-1)){
        if(length(sdG)>idNumThr){
            cat(paste("After filtering based on parameters 'naPer', 'meanPer' and 'varPer', the number of remaining ids is still larger than parameter 'idNumThr'. Thus, the function only selects the top ",idNumThr," ids with the largest variance\n",sep=""))
            sdG <- sdG[c(1:idNumThr)]
        }
    }
    
    cat("After removing ids based on parameters 'naPer', 'meanPer' and 'varPer', ",length(sdG)," ids are remained!\n",sep="")
    
    
    filterData <- filterData[names(sdG),]
    rm(meanG,sdG,inputMat)
    gc()
    
    enableWGCNAThreads()
	
    
	#calculate  correlation
	cat("Calculating ",corrType," correlation for each pair of ids...\n",sep="")
	threNum <- WGCNAnThreads()
	
	if(threNum>nThreads){
		threNum <- nThreads
	}
    
    if(missing==0){
        corrM <- WGCNA::cor(t(filterData),nThreads=threNum,method=corrType)
    }else{
        corrM <- WGCNA::cor(t(filterData),nThreads=threNum,use="pairwise.complete.obs",method=corrType)
    }
    
	disableWGCNAThreads()
    
    corrM[is.na(corrM)] <- 0
    corrM[corrM==1] <- 0.99999999999999999999
    corrM[corrM==(-1)] <- (-0.99999999999999999999)
    
    rankBest <- round(rankBest*nrow(corrM))
    if(rankBest==0){
        rankBest <- 1
    }
    
    if(matNetMethod=="value"){
        sigNetwork <- .valueBasedCoexp(corrM,valueThr,networkType)
        cat("Base on value-based method, an undirected network with ",length(union(sigNetwork[,1],sigNetwork[,2]))," nodes and ",nrow(sigNetwork)," edges was identified under threshold ",valueThr,".\n",sep="")
    }
    
    if(matNetMethod=="rank"){
        sigNetwork <- .rankBasedCoexp(filterData,corrM,rankBest,matNetMethod,networkType,netFDRMethod,netFDRThr)
        cat("Base on rank-based method, an undirected network with ",length(union(sigNetwork[,1],sigNetwork[,2]))," nodes and ",nrow(sigNetwork)," edges was identified when selecting most similar ",rankBest," nodes.\n",sep="")
    }
    
    if(matNetMethod=="directed"){
        sigNetwork <- .rankBasedCoexp(filterData,corrM,rankBest=1,matNetMethod,networkType,netFDRMethod,netFDRThr)
        cat("Base on directed-based method, a directed network with ",length(union(sigNetwork[,1],sigNetwork[,2]))," nodes and ",nrow(sigNetwork)," edges was identified when only selecting most similar node. \n",sep="")

    }
    
    return(sigNetwork)
}


.valueBasedCoexp <- function(corrArray,valueThr,networkType){
    
    diag(corrArray) <- 0
    
    if(networkType=="unsigned"){
        corrArray <- abs(corrArray)
    }
    
    corrArray[lower.tri(corrArray)] <- 0
    
    sigIndex <- which(corrArray>valueThr,arr.ind=TRUE)
    
    sigPair <- array(rownames(corrArray)[sigIndex],dim=dim(sigIndex))
    return(sigPair)
}

.rankBasedCoexp <- function(filterData,corrArray,rankBest,matNetMethod,networkType,netFDRMethod,netFDRThr){
    
    
    corrArray <- .calculateSigCorr(filterData,corrArray,netFDRMethod,netFDRThr)
    
    
    #calculatePosition <- function(sindex,rownum){
    #     col <- ceiling(sindex/rownum)
    #       correctI <- sindex-cumsum(1:col)[col]
    #       return(correctI)
    #   }
    
    
    if(networkType=="unsigned"){
        corrArray <- abs(corrArray)
    }
    
    diag(corrArray) <- 0
    
    bestNeighbor <- apply(corrArray,1,.findBestNeighbor,rankBest)
    
    corrArray[corrArray!=0] <- 0
    
    
    for(i in c(1:nrow(corrArray))){
        if(!is.null(bestNeighbor[[i]][[1]])){
            x <- bestNeighbor[[i]][[1]]
            corrArray[i,x] <- 1
        }
    }
    
    if(matNetMethod=="rank"){
        corrArray <- corrArray+t(corrArray)
        corrArray[lower.tri(corrArray)] <- 0
        sigIndex <- which(corrArray==2,arr.ind=TRUE)
    }
    
    if(matNetMethod=="directed"){
        sigIndex <- which(corrArray==1,arr.ind=TRUE)
    }
    
    sigPair <- array(rownames(corrArray)[sigIndex],dim=dim(sigIndex))

    return(sigPair)
}

.calculateSigCorr <- function(filterData,corrArray,netFDRMethod,netFDRThr){
    n <- t(!is.na(t(filterData))) %*% (!is.na(t(filterData)))
    suppressWarnings(t <- (corrArray * sqrt(n - 2))/sqrt(1 - corrArray^2))
    suppressWarnings(corrP <- 2 * (1 - pt(abs(t), (n - 2))))
    rm(n,t)
    gc()
    corrP[is.na(corrP)] <- 0
    corrP[corrP>1] <- 1
    diag(corrP) <- (-1)
    corrSig <- apply(corrP,1,.fdrFilter,netFDRMethod,netFDRThr)
    corrSig <- t(corrSig)
    rm(corrP)
    gc()
    corrArray <- corrArray*corrSig
    return(corrArray)

}

.fdrFilter <- function(vector,netFDRMethod,netFDRThr){
    a <- vector[vector>=0]
    adjp <- p.adjust(a,method=netFDRMethod)
    vector[names(adjp)] <- adjp
    vector[vector==(-1)] <- 1
    sigP <- ifelse(vector<=netFDRThr,1,0)
    return(sigP)
}



.findBestNeighbor <- function(corrList,rankBest){
    
    corrList <- corrList[corrList!=0]
    
    if(length(corrList)==0){
        return(NULL)
    }else{
    
        if(length(corrList)<rankBest){
            return(list(names(corrList)))
        }else{
            corrList <- sort(corrList,decreasing=T)
            return(list(names(corrList)[1:rankBest]))
        }
    }
}

.calculateNANum <- function(vector){
    return(sum(is.na(vector))/length(vector))
}