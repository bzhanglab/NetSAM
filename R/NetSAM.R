NetSAM <-
function(inputNetwork, outputFileName, outputFormat="nsm", edgeType="unweighted", map_to_genesymbol=FALSE, organism="hsapiens", idType="auto", minModule=0.003, stepIte=FALSE, maxStep=4, moduleSigMethod="cutoff", modularityThr=0.2, ZRanNum=10, PerRanNum=100, ranSig=0.05, edgeThr=(-1), nodeThr=(-1), nThreads=3){
	
    
    #require(igraph) || stop("Package igraph version 0.6 is required!")
	#require(seriation) || stop("Package seriation version 1.0-10 is required!")
	#require(graph) || stop("Package graph version 1.40.1 is required!")
    #require(WGCNA) || stop("Package WGCNA version 1.34.0 is required!")
    #require(doSNOW) || stop("Package doSNOW version 1.0.6 is required!")
    #require(foreach) || stop("Package foreach version 1.4.0 is required!")
    #require(multicore) || stop("Package multicore version 0.1-7 is required!")
    #require(tools) || stop("Package tools version 3.0.0 is required!")
    
    organisms <- c("hsapiens","mmusculus","rnorvegicus","drerio","celegans","scerevisiae","cfamiliaris","dmelanogaster","athaliana")
    names(organisms) <- c("Hs","Mm","Rn","Dr","Ce","Sc","Cf","Dm","At")
    
    if(length(which(organisms==organism))==0){
        stop("Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana!")
    }
    
    
	if(missing(inputNetwork)){
		stop("Please input the network!\n")
	}
    
    
    #load Network
	if(class(inputNetwork)=="character"){
        if(file_ext(inputNetwork)!="net"){
            stop("The extension of the input file should be 'net'!\n")
        }else{
            network <- read.graph(inputNetwork,format="ncol")
            if(edgeType=="weighted"){
                if(is.null(E(network)$weight)){
                    stop("The input network does not contain edge weights. Please add the edge weights or change parameter 'edgeType' to 'unweigthed'!")
                }
            }else{
                if(!is.null(E(network)$weight)){
                    stop("The input network contains edge weights. Please remove the edge weights or change parameter 'edgeType' to 'weigthed'!")
                }
            }
        }
	}else{
		if(class(inputNetwork)=="data.frame" || class(inputNetwork)=="matrix"){
            if(edgeType=="weighted"){
                if(ncol(inputNetwork)!=3){
                    stop("Data object should contain three columns: interactor1, interactor2 and edge weight!")
                }else{
                    weight <- inputNetwork[,3]
                    if(!is.numeric(weight)){
                        stop("The edge weight should be numeric!")
                    }else{
                        inputNetwork <- as.matrix(inputNetwork[,c(1,2)])
                        inputNetwork_S <- as.character(inputNetwork)
                        inputNetwork_S <- array(inputNetwork_S,dim=dim(inputNetwork))
                        network <- graph.edgelist(inputNetwork_S,directed=F)
                        E(network)$weight <- weight
                    }
                }
            }else{
                if(ncol(inputNetwork)!=2){
                    stop("data object should contain two columns!\n");
                }else{
                    inputNetwork <- as.matrix(inputNetwork)
                    inputNetwork_S <- as.character(inputNetwork)
                    inputNetwork_S <- array(inputNetwork_S,dim=dim(inputNetwork))
                    network <- graph.edgelist(inputNetwork_S,directed=F)
                }
            }
            rm(inputNetwork,inputNetwork_S)
            gc()
		}else{
			if(class(inputNetwork)=="igraph"){
                if(edgeType=="unweighted"){
                    if(!is.null(E(inputNetwork)$weight)){
                        stop("The input network contains edge weights. Please remove the edge weights or change parameter 'edgeType' to 'weigthed'!")
                    }else{
                        network <- inputNetwork
                    }
                }else{
                    if(is.null(E(inputNetwork)$weight)){
                        stop("The input network does not contain edge weights. Please add the edge weights or change parameter 'edgeType' to 'unweigthed'!")
                    }else{
                        network <- inputNetwork
                    }
                }
                rm(inputNetwork)
                gc()
            }else{
                if(class(inputNetwork)=="graphNEL"){
                    network <- igraph.from.graphNEL(inputNetwork)
                    if(edgeType=="unweighted"){
                        if(!is.null(E(network)$weight)){
                            stop("The input network contains edge weights. Please remove the edge weights or change parameter 'edgeType' to 'weigthed'!")
                        }
                    }else{
                        if(is.null(E(network)$weight)){
                            stop("The input network does not contain edge weights. Please add the edge weights or change parameter 'edgeType' to 'unweigthed'!")
                        }
                    }
					rm(inputNetwork)
					gc()
				}else{
                    stop("The input network should be from a file or a data object with data.frame, matrix, graphNEL or igraph class. Other types of input are invalid!\n")
				}

            }
        }
	}
	
	if(outputFormat != "none" && missing(outputFileName)){
		stop("Please input the output file name!\n")
	}else{
		if(substr(outputFileName,nchar(outputFileName),nchar(outputFileName))=="/"){
			stop("Please input the output file name!\n")
		}
	}
	
    if(length(which(outputFormat %in% c("nsm","gmt","multiple","none")))==0){
		stop("The input 'outputFormat' is invalid! Please select an output format from 'nsm', 'gmt', 'multiple' and 'none'!\n")
	}
    
    if(!is.logical(map_to_genesymbol)){
        stop("The input 'map_to_genesymbol' is invalid! 'map_to_genesymbol' should only be TRUE or FALSE!\n")
    }
    
    if(length(which(organisms==organism))==0){
        stop("Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana!")
    }
    
    if(minModule<0 || minModule>0.2){
        stop("The input 'minModule' is invalid! 'minModule' should be from 0 to 0.2!\n")
    }
    

    if(maxStep<2 || class(maxStep)!="numeric"){
        stop("The input 'maxStep' is invalid! 'maxStep' should be a number and larger than 2!\n")
    }
    

	if(length(which(moduleSigMethod %in% c("cutoff","zscore","permutation")))==0){
		stop("The input 'moduleSigMethod' is invalid! Please select an option from 'cutoff','zscore' and 'permutation'!\n")
	}
    
    if(modularityThr>1 || modularityThr<0 || class(modularityThr)!="numeric"){
        stop("The input 'modularityThr' is invalid! 'modularityThr' should be a positive number and less than 1!\n")
    }
    
    if(ZRanNum<10 || class(ZRanNum)!="numeric"){
        stop("The input 'ZRanNum' is invalid! 'ZRanNum' should be a number and larger than 10!\n")
    }
    
    
    if(PerRanNum<10 || class(PerRanNum)!="numeric"){
        stop("The input 'PerRanNum' is invalid! 'PerRanNum' should be a number and larger than 10!\n")
    }
    
    
    if(ranSig>1 || ranSig <0 || class(ranSig)!="numeric"){
        stop("The input 'ranSig' is invalid! 'ranSig' should be a positive number and less than 1!\n")
    }
    
    if(edgeThr<(-1) || class(edgeThr)!="numeric"){
        stop("The input 'edgeThr' is invalid! 'edgeThr' should be -1 or a positive number!\n")
    }
    
    if(nodeThr<(-1) || class(nodeThr)!="numeric"){
        stop("The input 'nodeThr' is invalid! 'nodeThr' should be -1 or a positive number!\n")
    }
    
    if(nThreads<1 || class(nThreads)!="numeric"){
        stop("The input 'nThreads' is invalid! 'nThreads' should be a number and greater than 1!\n")
    }
    
    
    if(!is.simple(network)){
        warning("The input network contain loop or multiple edges. For unweighted network, NetSAM will remove loop or multiple edges. For weighted network, NetSAM will sum the weights from loop or multiple edges!!!")
        network <- simplify(network)
    }
    
    if(nodeThr!=(-1)){
        if(vcount(network)>nodeThr){
            stop("The number of nodes in the network is over the threshold ",nodeThr,". Please adjust the parameters to decrease the network size!\n")
        }
    }
    
    if(edgeThr!=(-1)){
        if(ecount(network)>edgeThr){
            stop("The number of edges in the network is over the threshold ",edgeThr,". Please adjust the parameters to decrease the network size!\n")
        }
    }


    ###################Get the whole gene symbols for the given organism###################
    if(organism=="scerevisiae"){
        #require("org.Sc.sgd.db") || stop("package org.Sc.sgd.db is required!!\n")
        genesymbol <- keys(eval(parse(text="org.Sc.sgd.db::org.Sc.sgdCOMMON2ORF")))
    }else{
        if(organism=="athaliana"){
            #require("org.At.tair.db") || stop("package org.At.tair.db is required!!\n")
            genesymbol <- unlist(as.list(eval(parse(text="org.At.tair.db::org.At.tairSYMBOL"))))
            names(genesymbol) <- NULL
            genesymbol <- unique(genesymbol)
            genesymbol <- genesymbol[!is.na(genesymbol)]
            
        }else{
            or <- names(organisms)[which(organisms==organism)]
            #require(paste("org.",or,".eg.db",sep=""),character.only = TRUE) || stop(paste("package org.", or, ".eg.db is required!!\n", sep=""))
            genesymbol <- keys(eval(parse(text=paste("org.",or,".eg.db::","org.",or,".egSYMBOL2EG",sep=""))))
        }
    }
    ########################################################################################
    
    geneid <- V(network)$name
    
    is.genesymbol <- FALSE
    if(length(intersect(geneid,genesymbol))>0.6*length(geneid)){
        is.genesymbol <- TRUE
    }else{
        is.genesymbol <- FALSE
    }
    
    mapping.status <- 0

    if(map_to_genesymbol==TRUE){
        ne <- mapToSymbol(inputData=network,organism=organism,inputType="network",idType=idType, edgeType=edgeType, verbose=TRUE)
        if(!is.null(ne)){
            network_edgelist <- ne$network
            netMat <- as.matrix(network_edgelist[,c(1,2)])
            network <- graph.edgelist(netMat,directed=FALSE)
            if(edgeType=="weighted"){
                E(network)$weight <- network_edgelist[,3]
            }
            mapping.status <- 0
        }else{
            stop("Please check the parameter 'idType' or set the parameter 'map_to_genesymbol' as FALSE!\n")
        }
    }else{
        if(is.genesymbol==TRUE){
            mapping.status <- 0
        }else{
            if(idType==""){     #NONE OF THE ABOVE
                mapping.status <- 1
            }else{
                ne <- mapToSymbol(inputData=network,organism=organism,inputType="network",idType=idType,edgeType=edgeType,verbose=FALSE)
                if(!is.null(ne)){
                    idmap <- ne$idmap
                    mapping.status <- 2
                }else{
                    mapping.status <- 1
                }
            }
        }
    }
    
    minModule <- round(vcount(network)*minModule)
    if(minModule<5){
        minModule <- 5
    }
    
    enableWGCNAThreads()
    threNum <- WGCNAnThreads()
    if(threNum>nThreads){
        threNum <- nThreads
    }
    


    network_module <- .identifyModule(network,minModule,stepIte,maxStep,moduleSigMethod,modularityThr,ZRanNum,PerRanNum,ranSig,threNum,edgeType)
    
    geneorder <- network_module$geneorder
    
    disableWGCNAThreads()
    
	if(outputFormat=="nsm" || outputFormat=="none" || outputFormat=="multiple"){
        
        network_edges <- get.edgelist(network)
        if(edgeType=="weighted"){
            if(ecount(network)>500000){
                warning("The input weighted network contain over 500,000 edges. Because nsm file contains edge information, remaining all edges of the network will cause a very large size of nsm file, which will be difficult to upload to NetGestalt. NetSAM only remains top 500,000 edges based on the edge weights from largest to smallest!!\n")
                w <- E(network)$weight
                w <- order(w,decreasing=T)
                network_edges <- network_edges[w,]
                network_edges <- network_edges[1:500000,]
            }
        }
        
        hmiFile <- .createHMIFile(network_module$allHir,geneorder)
        
        rownames(geneorder) <- geneorder[,4]
        netgestalt <- list(rulfile=geneorder,hmifile=hmiFile,network=network_edges)
        
        if(outputFormat=="multiple"){
            outputFile <- paste(outputFileName,"_rulerfile.rul",sep="")
            write.table(geneorder,file=outputFile,row.names=F,col.names=T,sep="\t",quote=F)
            
            outputFile <- paste(outputFileName,"_hmifile.hmi",sep="")
            write.table(hmiFile,file=outputFile,row.names=F,col.names=T,sep="\t",quote=F)
            
            outputFile <- paste(outputFileName,"_network.net",sep="")
            write.table(network_edges,file=outputFile,row.names=F,col.names=F,sep="\t",quote=F)
            
        }
        
        if(outputFormat=="nsm"){
            outputFile <- paste(outputFileName,".nsm",sep="")
            note <- "##  Ruler file  ##"
            write.table(note,file=outputFile,row.names=F,col.names=F,quote=F)
            suppressWarnings(write.table(geneorder,file=outputFile,row.names=F,col.names=T,append=T,quote=F,sep="\t"))
            note <- "##  HMI file ##"
            write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
            suppressWarnings(write.table(hmiFile,file=outputFile,row.names=F,col.names=T,quote=F,append=T,sep="\t"))
            note <- "##  Network file  ##"
            write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
            write.table(network_edges,file=outputFile,row.names=F,col.names=F,quote=F,append=T,sep="\t")
            
            netgestalt$mapping.status <- mapping.status
            
            note <- "##  Mapping Status  ##"
            write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
            suppressWarnings(write.table(paste("mapping.status=",mapping.status,sep=""),file=outputFile,row.names=F,col.names=F,quote=F,append=T,sep="\t"))
            
            
            if(mapping.status==2){
                note <- "##  Gene Annotation  ##"
                write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
                suppressWarnings(write.table(idmap,file=outputFile,row.names=F,col.names=T,quote=F,append=T,sep="\t"))
                netgestalt$geneAnn <- idmap
            }

		}
        
        cat("NetSAM identified ",(nrow(hmiFile)-1)," modules in ",(length(unique(hmiFile[,2]))-1)," levels!\n",sep="")
        cat("Processing completed!\n\n")
        return(netgestalt)
        
        
    }else{
        hmiFile <- .createHMIFile(network_module$allHir,geneorder)
        gmtFile <- .createGMTFile(hmiFile,geneorder)
        outputFile1 <- paste(outputFileName,".gmt",sep="")
        write.table(gmtFile,file=outputFile1,row.names=F,col.names=F,sep="\t",quote=F)
        result <- gmtFile
        cat("NetSAM identified ",(nrow(hmiFile)-1)," modules in ",(length(unique(hmiFile[,2]))-1)," levels!\n",sep="")
        cat("Processing completed!\n\n")
        return(result)
    }
	
}


.identifyModule <-
function(network, minModule, stepIte, maxStep, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum, edgeType){
	
    cat("\nIdentifying the hierarchical modules of the network...\n")
    overlap_genesymbol_networkP <- V(network)$name
    subnetworkInfo <- list()
    subnetworkInfo_index <- 1
    network_cluster <- clusters(network)
    network_cluster_size <- network_cluster$csize
        
    if(length(network_cluster_size[network_cluster_size>=minModule])==0){
        stop("The size of all subnetworks are less than ",minModule,". Please adjust the parameter 'minModule'!\n\n")
    }
        
    network_cluster_size <- data.frame(id=c(1:length(network_cluster_size)),cluster_size=network_cluster_size,stringsAsFactors=F)
    network_cluster_size <- network_cluster_size[order(-network_cluster_size[,2]),]
    network_cluster_membership <- network_cluster$membership
        
    network_ordered_node <- vector()
    subnetwork_id <- 1
    for(i in c(1:nrow(network_cluster_size))){
        sub_network_size <- network_cluster_size[i,2]
        if(sub_network_size>=minModule){
            cat("Starting to analysis connected component ",subnetwork_id,"!\n",sep="")
            subnetwork_id <- subnetwork_id+1
            subnetwork_node <- overlap_genesymbol_networkP[which(network_cluster_membership==network_cluster_size[i,1])]
            subnetwork <- .extractSubNetwork(network,subnetwork_node,edgeType)
            allHir <- .identifyHierOr(subnetwork, minModule, stepIte, maxStep, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum, edgeType)
            subnetworkInfo[[subnetworkInfo_index]] <- allHir
            subnetworkInfo_index <- subnetworkInfo_index+1
            network_ordered_node <- c(network_ordered_node,subnetwork_node)
            cat("\n")
        }
            
    }
        
        
    cat("\nReordering the genes in the one dimentional layout...\n")
        
    geneorder <- data.frame(ruler_id=c(1:length(overlap_genesymbol_networkP)),node_type="Gene",node_db="Entrez Gene",node_db_id=overlap_genesymbol_networkP,node_name=overlap_genesymbol_networkP,stringsAsFactors=F)
    rownames(geneorder) <- overlap_genesymbol_networkP
    unann_node <- setdiff(overlap_genesymbol_networkP,network_ordered_node)
    if(length(unann_node) !=0){
        network_ordered_node <- c(network_ordered_node,unann_node)
    }
    geneorder <- geneorder[network_ordered_node,]
        
        
    geneorder <- .orderDiffLevel(subnetworkInfo,geneorder)
        
    result <- list(network=network,allHir=subnetworkInfo,geneorder=geneorder)
    return(result)
}


.calculateRandomWalkerAdjectMatrix <- function(network_igraph,step,edgeType){
    
	#calculate the random walk distance for the network
    protein_in_ppi <- V(network_igraph)$name
    if(edgeType=="unweighted"){
        E(network_igraph)$weight <- rep(1,ecount(network_igraph))
    }
    
    nweight <- graph.strength(network_igraph,mode="all")
    V(network_igraph)$weight <- as.numeric(nweight)
    
    
    eweight <- E(network_igraph)$weight
    
    degree_node <- igraph::degree(network_igraph)
    network_igraph <- add.edges(network_igraph,rbind(protein_in_ppi,protein_in_ppi)) #walktrap algorithm adds all self interactions for each node
    E(network_igraph)$weight <- c(eweight,V(network_igraph)$weight/degree_node)
    V(network_igraph)$weight <- V(network_igraph)$weight+ V(network_igraph)$weight/degree_node
    
    nweight <- V(network_igraph)$weight
    
    adjMatrix <- get.adjacency(network_igraph,attr="weight")
    adjMatrix <- as.matrix(adjMatrix)
    
    W_weighted <- adjMatrix/nweight
    
    rm(adjMatrix,network_igraph,eweight,degree_node)
    gc()
    
    #cat("Create transition Matrix...\n")
    
    tranM <- W_weighted
    
    if(step>1){
        for(i in c(1:(step-1))){
            tranM <- tranM%*%W_weighted
        }
    }
    
    rm(W_weighted)
    gc()
    
    #tranM <- as.matrix(tranM)
    
    tranM <- t(tranM)/sqrt(nweight)
    tranM <- t(tranM)
    
    #cat("Calculate EU distance...\n")
    
    smat <- apply(tranM,1,crossprod)
    mat1 <- matrix(smat,nrow=length(protein_in_ppi),ncol=length(protein_in_ppi))
    mat3 <- tcrossprod(tranM)
    mat4 <- mat1+t(mat1)-2*mat3
    rm(smat,mat1,mat3,tranM)
    gc()
    mat4[mat4<0] <- 0
    diag(mat4) <- 0
    mat4 <- sqrt(mat4)
    
    return(mat4)
}


.transformFromWalktrapToHclust <- function(walktrap){
    
    #transform walktrap to hclust
    
    merge <- walktrap$merges
    name <- walktrap$names
    N <- length(name)
    
    merge[which(merge[,1]<=N),1] <- (-merge[which(merge[,1]<=N),1])
    merge[which(merge[,2]<=N),2] <- (-merge[which(merge[,2]<=N),2])
    
    merge[which(merge[,1]>N),1] <- (merge[which(merge[,1]>N),1]-N)
    merge[which(merge[,2]>N),2] <- (merge[which(merge[,2]>N),2]-N)
    
    height <- c(1:(N-1))
    
    dend <- as.dendrogram(walktrap)
    order <- order.dendrogram(dend)
    
    labels <- name
    
    ht <- list(merge=merge,height=height,order=order,labels=labels,method="walktrap",call=match.call(),dist.method="randomwalk")
    
    class(ht) <- "hclust"
    
    return(ht)
}

.walktrapcommunity <- function(network_igraph,steps,edgeType){
    if(edgeType=="weighted"){
        network_walktrap <- walktrap.community(network_igraph,weights=E(network_igraph)$weight,steps=steps)
    }else{
        network_walktrap <- walktrap.community(network_igraph,steps=steps)
    }
    return(network_walktrap)
}

.evaluateWalktrapStep <- function(network_igraph,stepIte,maxStep,level,threNum,edgeType,subM){
    
	#evaluate the optimal Step for the network
    
    if(stepIte==FALSE || maxStep==2 || level==1){
        network_walktrap <- .walktrapcommunity(network_igraph,steps=maxStep,edgeType)
        network_adjMatrix <- .calculateRandomWalkerAdjectMatrix(network_igraph,maxStep,edgeType)
        rm(network_igraph)
        gc()
        network_adjMatrix <- as.dist(network_adjMatrix)
        network_hclust <- .transformFromWalktrapToHclust(network_walktrap)
        network_order <- seriate(network_adjMatrix,method="OLO",control=list(hclust=network_hclust))
        network_order <- get_order(network_order)
        network_info <- list(gene=network_walktrap$names,membership=network_walktrap$membership,modularity=max(network_walktrap$modularity),step=maxStep,order=network_order,level=level)
        return(network_info)
    }else{
        
        if(maxStep < threNum){
            threNum <- maxStep
            perN <- 1
        }else{
            perN <- floor((maxStep-1)/threNum)
        }

        
        cl <- makeCluster(threNum)
        registerDoSNOW(cl)

        
        i <- 1
        optimalM <- foreach(i=1:threNum, .packages="igraph") %dopar% {
            finalStep <- data.frame(step=0,modularity=1,stringsAsFactors=F)
            fi <- 1
            if(i==1){
                step_start <- 2
                step_end <- i*perN+1
            }else{
                step_start <- i*perN-(perN-2)
                if(i==threNum){
                    step_end <- maxStep
                }else{
                    step_end <- i*perN+1
                }
            }
            for(j in c(step_start:step_end)){
                if(edgeType=="weighted"){
                    w <- walktrap.community(network_igraph,weights=E(network_igraph)$weight,steps=j)
                }else{
                    w <- walktrap.community(network_igraph,steps=j)
                }
                m <- max(w$modularity)
                finalStep[fi,1] <- j
                finalStep[fi,2] <- m
                fi <- fi + 1
            }
            
            return(finalStep)
        }
        stopCluster(cl)
        
        finalStep <- data.frame(step=0,modularity=1,stringsAsFactors=F)
        for(i in c(1:length(optimalM))){
            finalStep <- rbind(finalStep,optimalM[[i]])
        }
        finalStep <- finalStep[finalStep[,1]!=0,]
        finalStep <- finalStep[order(-finalStep[,2],finalStep[,1]),]
        
        optimalStep <- finalStep[1,1]
        cat("The optimal Step of network ",subM," is ",optimalStep," !\n",sep="")
        
        optimalwalktrap <- .walktrapcommunity(network_igraph,steps=optimalStep,edgeType)
        network_adjMatrix <- .calculateRandomWalkerAdjectMatrix(network_igraph,optimalStep,edgeType)
        rm(network_igraph)
        gc()
        network_adjMatrix <- as.dist(network_adjMatrix)
        network_hclust <- .transformFromWalktrapToHclust(optimalwalktrap)
        network_order <- seriate(network_adjMatrix,method="OLO",control=list(hclust=network_hclust))
        network_order <- get_order(network_order)
        maxWalktrap <- list(gene=optimalwalktrap$names,membership=optimalwalktrap$membership,modularity=max(optimalwalktrap$modularity),step=optimalStep,order=network_order,level=level)
        return(maxWalktrap)
    }
}


.identifySig_Unweighted <- function(network_igraph,network_modularity,step,moduleSigMethod,modularityThr, ZRanNum, PerRanNum, ranSig, threNum){
    
	#identify whether the network can be separated again
    
    degree <- igraph::degree(network_igraph)
    sig <- 0
    i <- 1
    
    if(moduleSigMethod=="cutoff"){
        if(network_modularity>modularityThr){
            sig <- 1
        }
    }
    
    if(moduleSigMethod=="zscore"){
        
        if(ZRanNum<threNum){
            threNum <- ZRanNum
            perN <- 1
        }else{
            perN <- floor(ZRanNum/threNum)
        }


        cl <- makeCluster(threNum)
        registerDoSNOW(cl)
        
        ranmodu <- foreach(i=1:threNum, .combine="c", .packages="igraph") %dopar%{
            
            randomM <- c()
            
            rN <- perN
            if(i==threNum){
                rN <- ZRanNum-(i-1)*perN
            }
            
            for(j in c(1:rN)){
                suppressWarnings(rannet <- degree.sequence.game(degree,method="vl"))
                ran_walktrap <- walktrap.community(rannet,steps=step)
                ranModularity <- max(ran_walktrap$modularity)
                randomM <- c(randomM,ranModularity)
            }
            return(randomM)
        }
        stopCluster(cl)
        
        ranmodu_mean <- mean(ranmodu)
        ranmodu_sd <- sd(ranmodu)
        if(ranmodu_sd == 0){
            if(network_modularity>ranmodu_mean){
                sig <- 1
            }
        }else{
            Z_network_modularity <- (network_modularity-ranmodu_mean)/ranmodu_sd
            p <- pnorm(Z_network_modularity,mean=0,sd=1,lower.tail=F)
            if(p<ranSig){
                sig <- 1
            }
        }
    }
    
    if(moduleSigMethod=="permutation"){
        
        if(PerRanNum<threNum){
            threNum <- PerRanNum
            perN <- 1
        }else{
            perN <- floor(PerRanNum/threNum)
        }

        
        cl <- makeCluster(threNum)
        registerDoSNOW(cl)
        
        ranmodu <- foreach(i=1:threNum, .combine="c", .packages="igraph") %dopar%{
            
            randomM <- c()
            
            rN <- perN
            if(i==threNum){
                rN <- PerRanNum-(i-1)*perN
            }
            
            for(j in c(1:rN)){
            
                suppressWarnings(rannet <- degree.sequence.game(degree,method="vl"))
                ran_walktrap <- walktrap.community(rannet,steps=step)
                ranModularity <- max(ran_walktrap$modularity)
                randomM <- c(randomM,ranModularity)
            }
            return(randomM)
        }
        stopCluster(cl)
        
        p <- length(ranmodu[ranmodu>=network_modularity])/PerRanNum
        if(p<ranSig){
            sig <- 1
        }
    }
    
    return(sig)
}


.identifySig_Weighted <- function(network_igraph, network_modularity, step, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum){
    
	#identify whether the network can be separated again
    
    
    sig <- 0
    i <- 1
    if(moduleSigMethod=="cutoff"){
        if(network_modularity>modularityThr){
            sig <- 1
        }
    }
    
    if(moduleSigMethod=="zscore"){
        
        #network_igraph <- network_info$network
        weight <- E(network_igraph)$weight
        #ranmodu <- vector()
        #step <- network_info$step
        
        if(ZRanNum<threNum){
            threNum <- ZRanNum
            perN <- 1
        }else{
            perN <- floor(ZRanNum/threNum)
        }
        
        cl <- makeCluster(threNum)
        registerDoSNOW(cl)
        
        ranmodu <- foreach(i=1:threNum, .combine="c", .packages="igraph") %dopar%{
            
            randomM <- c()
            
            rN <- perN
            if(i==threNum){
                rN <- ZRanNum-(i-1)*perN
            }
            
            for(j in c(1:rN)){
            
                ranWeight <- sample(weight,length(weight))
                ranModularity <- max(walktrap.community(network_igraph,weights=ranWeight,steps=step)$modularity)
                randomM <- c(randomM,ranModularity)
            }
            return(randomM)
        }
        stopCluster(cl)
        rm(network_igraph,weight)
        gc()
        #cat("randomMo:",ranmodu,"\n")
        ranmodu_mean <- mean(ranmodu)
        ranmodu_sd <- sd(ranmodu)
        if(ranmodu_sd == 0){
            if(network_modularity>ranmodu_mean){
                sig <- 1
            }
        }else{
            Z_network_modularity <- (network_modularity-ranmodu_mean)/ranmodu_sd
            p <- pnorm(Z_network_modularity,mean=0,sd=1,lower.tail=F)
            if(p<ranSig){
                sig <- 1
            }
        }
    }
    
    if(moduleSigMethod=="permutation"){
        
        #network_igraph <- network_info$network
        weight <- E(network_igraph)$weight
        #step <- network_info$step
        
        if(PerRanNum<threNum){
            threNum <- PerRanNum
            perN <- 1
        }else{
            perN <- floor(PerRanNum/threNum)
        }
       
        cl <- makeCluster(threNum)
        registerDoSNOW(cl)
        
        ranmodu <- foreach(i=1:threNum, .combine="c", .packages="igraph") %dopar%{
            
            randomM <- c()
            
            rN <- perN
            if(i==threNum){
                rN <- PerRanNum-(i-1)*perN
            }
            
            for(j in c(1:rN)){
            
                ranWeight <- sample(weight,length(weight))
                ranModularity <- max(walktrap.community(network_igraph,weights=ranWeight,steps=step)$modularity)
                randomM <- c(randomM,ranModularity)
            }
            return(randomM)
        }
        stopCluster(cl)
        rm(network_igraph,weight)
        gc()
        
        p <- length(ranmodu[ranmodu>=network_modularity])/PerRanNum
        if(p<ranSig){
            sig <- 1
        }
    }
    
    return(sig)
}


.identifySig <- function(network_igraph, network_modularity, step, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum, edgeType){
    if(edgeType=="weighted"){
        sig <- .identifySig_Weighted(network_igraph, network_modularity, step, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum)
    }else{
        sig <- .identifySig_Unweighted(network_igraph, network_modularity, step, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum)
    }
    
    return(sig)
}


.identifyHierOr <- function(network_maxComponent, minModule, stepIte, maxStep, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum, edgeType){
    
	#identify hierarchical organization of network
    
    allHir <- list()
    ai <- 1
    sigHir <- vector()
    si <- 1
    
    start <- 1
    
    cat("Evaluating networks in Level 1 ...\n")
    
    levelid <- 2
    
    network_info <- .evaluateWalktrapStep(network_maxComponent, stepIte, maxStep, level=1, threNum, edgeType, subM=1)
    cat("Network modularity: ",network_info$modularity,"\n\n",sep="")
    
    network_sig <- .identifySig(network_maxComponent, network_info$modularity, network_info$step, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum, edgeType)
    
    allHir[[ai]] <- network_info
    rm(network_info)
    gc()
    
    
    if(network_sig==1){
        sigHir <- c(sigHir,ai)
        si <- si+1
    }
    
    ai <- ai+1
    
    while(length(sigHir)>0 && start!=si){
        currentPos <- sigHir[start]
        currentNetwork <- allHir[[currentPos]]
        start <- start+1
        
        currentNetwork_level <- currentNetwork$level
        currentNetwork_node <- currentNetwork$gene
        
        currentNetwork_membership <- currentNetwork$membership
        currentNetwork_membership <- data.frame(node=currentNetwork_node,membership=currentNetwork_membership,stringsAsFactors=F)
        currentNetwork_membership_count <- tapply(currentNetwork_membership[,1],currentNetwork_membership[,2],length)
        currentNetwork_membership_count <- currentNetwork_membership_count[currentNetwork_membership_count>=minModule]
        
        if(length(currentNetwork_membership_count)==0){
            next
        }
        
        currentNetwork_membership_group <- as.integer(names(currentNetwork_membership_count))
        rm(currentNetwork_membership_count,currentNetwork_node)
        gc()
        
        subLevel <- currentNetwork_level+1
        if(subLevel==levelid){
            cat("Evaluating networks in Level ",levelid," ...\n",sep="")
            levelid <- levelid+1
            subM <- 1
        }
        for(i in c(1:length(currentNetwork_membership_group))){
            subnetwork_node <- currentNetwork_membership[currentNetwork_membership[,2]==currentNetwork_membership_group[i],1]
            subnetwork_igraph <- .extractSubNetwork(network_maxComponent,subnetwork_node,edgeType)
            subnetwork_info <- .evaluateWalktrapStep(subnetwork_igraph, stepIte, maxStep, subLevel, threNum, edgeType, subM)
            cat("Modularity of network ",subM,": ",subnetwork_info$modularity,"\n\n",sep="")
            subnetwork_sig <- .identifySig(subnetwork_igraph, subnetwork_info$modularity, subnetwork_info$step, moduleSigMethod, modularityThr, ZRanNum, PerRanNum, ranSig, threNum, edgeType)
            
            allHir[[ai]] <- subnetwork_info
            rm(subnetwork_node,subnetwork_igraph,subnetwork_info)
            gc()
            
            if(subnetwork_sig==1){
                sigHir <- c(sigHir,ai)
                si <- si+1
            }
            
            ai <- ai+1
            subM <- subM + 1
        }
    }
    
    return(allHir)
}


.orderDiffLevel <- function(subnetworkInfo,geneorder){
    
	#reorder all genes in the network according to the optimal position in each level
    
    for(l in c(1:length(subnetworkInfo))){
        allHir <- subnetworkInfo[[l]]
		
        for(i in c(1:length(allHir))){
            ori <- allHir[[i]]$order
            node <- allHir[[i]]$gene
            node <- node[ori]
            node <- data.frame(id=c(1:length(node)),name=node,stringsAsFactors=F)
            node <- node[order(node[,2]),]
            
            geneorder_subpos <- which(geneorder[,5] %in% node[,2])
            geneorder_sub <- geneorder[geneorder_subpos,]
            geneorder_sub <- geneorder_sub[order(geneorder_sub[,5]),]
            geneorder_sub <- cbind(geneorder_sub,node[,1])
            geneorder_sub <- geneorder_sub[order(geneorder_sub[,6]),]
            geneorder_sub <- geneorder_sub[,c(1:5)]
            geneorder[geneorder_subpos,] <- geneorder_sub
            
        }
		
    }
    geneorder[,1] <- c(1:nrow(geneorder))
    return(geneorder)
}

.createHMIFile <- function(subnetworkInfo,geneorder){
	#create HMI file
    hmiFile <- data.frame(best="N",level=0,order=1,name="ALL",start=1,end=nrow(geneorder),stringsAsFactors=F)
    hi <- 2
    for(l in c(1:length(subnetworkInfo))){
        allHir <- subnetworkInfo[[l]]
        for(i in c(1:length(allHir))){
            node <- allHir[[i]]$gene
            position <- which(geneorder[,5] %in% node)
            start <- min(position)
            end <- max(position)
            level <- allHir[[i]]$level
            if(level==2){
                best <- "Y"
            }else{
                best <- "N"
            }
            hmiFile[hi,1] <- best
            hmiFile[hi,2] <- level
            hmiFile[hi,3] <- 0
            hmiFile[hi,4] <- ""
            hmiFile[hi,5] <- start
            hmiFile[hi,6] <- end
            hi <- hi+1
        }
    }
    
    hmiFile <- hmiFile[order(hmiFile[,2],hmiFile[,5]),]
    allLevel <- sort(unique(hmiFile[,2]))
    for(i in c(1:length(allLevel))){
        position <- which(hmiFile[,2]==allLevel[i])
        start <- min(position)
        end <- max(position)
        hmiFile[start:end,3] <- c(1:(end-start+1))
    }
    hmiFile[,4] <- paste("Level",hmiFile[,2],"Module",hmiFile[,3],sep="_")
    return(hmiFile)
}


.createGMTFile <- function(hmiFile,geneorder){
    allgene <- geneorder[,4]
    allgene <- paste(allgene,collapse="\t")
    humancatfile <- data.frame(name="01",childnum=1,gene=allgene,start=1,end=nrow(geneorder),level=0,stringsAsFactors=F)
    
    for(i in c(2:nrow(hmiFile))){
        st <- hmiFile[i,5]
        en <- hmiFile[i,6]
        l <- hmiFile[i,2]
        for(j in c(1:nrow(humancatfile))){
            pl <- humancatfile[j,6]
            if(pl==(l-1)){
                ps <- humancatfile[j,4]
                pe <- humancatfile[j,5]
                if(st>=ps && en<=pe){
                    pcu <- humancatfile[j,1]
                    pcu_c <- humancatfile[j,2]
                    if(pcu_c < 10){
                        pcu_c <- paste("0",pcu_c,sep="")
                    }else{
                        pcu_c <- as.character(pcu_c)
                    }
                    
                    c_cu <- paste(pcu,pcu_c,sep="-")
                    humancatfile[j,2] <- humancatfile[j,2] + 1
                    humancatfile[i,1] <- c_cu
                    humancatfile[i,2] <- 1
                    cg <- geneorder[c(st:en),4]
                    cg <- paste(cg,collapse="\t")
                    humancatfile[i,3] <- cg
                    humancatfile[i,4] <- st
                    humancatfile[i,5] <- en
                    humancatfile[i,6] <- l
                    break
                }
            }
        }
    }
    humancatfile[,7] <- humancatfile[,5]-humancatfile[,4]+1
    gmtFile <- humancatfile[,c(1,7,3)]
    return(gmtFile)
}


.extractSubNetwork <- function(network,subList,edgeType){
    if(edgeType=="unweighted"){
        x <- get.edgelist(network)
        x <- x[x[,1] %in% subList & x[,2] %in% subList,]
        network <- graph.edgelist(x,directed=F)
    }else{
        x <- get.edgelist(network)
        weight <- get.edge.attribute(network,"weight")
        ind <- which(x[,1] %in% subList & x[,2] %in% subList)
        x <- x[ind,]
        weight <- weight[ind]
        network <- graph.edgelist(x,directed=F)
        E(network)$weight <- weight
    }
    return(network)
}




