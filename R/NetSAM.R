NetSAM <-
function(inputNetwork, outputFileName, minModule=(-1), maxStep=4, method="Modularity Cutoff", ModularityThr=0.2, ZRandomNum=10, permuteNum=100, pThr=0.05){
	
	calculateRandomWalkerAdjectMatrix <- function(network_igraph,step){
		
	#calculate the random walk distance for the network
		protein_in_ppi <- V(network_igraph)$name
		network_igraph <- add.edges(network_igraph,rbind(protein_in_ppi,protein_in_ppi)) #walktrap algorithm adds all self interactions for each node
		adjMatrix <- get.adjacency(network_igraph)
		degree <- igraph:::degree(network_igraph)-1    #after adding self interactions, the degree will increase 2 for each node. Thus, we should substract degree by 1.
		W_unweighted <- adjMatrix*(1/degree)
		
		#cat("Create transition Matrix...\n")
		
		tranM <- W_unweighted
		
		if(step>1){
			for(i in c(1:(step-1))){
				tranM <- tranM%*%W_unweighted
			}
		}
		
		rm(W_unweighted)
		
		degreesquar <- sqrt(degree)
		
		tranM <- as.matrix(tranM)
		
		tranMdegree <- t(tranM)/degreesquar
		tranMdegree <- t(tranMdegree)
		
		rm(tranM)
		gc()
		#cat("Calculate EU distance...\n")
		
		smat <- apply(tranMdegree,1,crossprod)
		mat1 <- matrix(smat,nrow=length(protein_in_ppi),ncol=length(protein_in_ppi))
		mat3 <- tcrossprod(tranMdegree)
		mat4 <- mat1+t(mat1)-2*mat3
		mat4[mat4<0] <- 0
		diag(mat4) <- 0
		adjMatrix <- sqrt(mat4)
		
		adjMatrix <- as.matrix(adjMatrix)
		
		return(adjMatrix)
	}
	
	transformFromWalktrapToHclust <- function(walktrap){
		
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
	
	
	evaluateWalktrapStep <- function(network_igraph,maxStep,level){
		
	#evaluate the optimal Step for the network
		
		network_info <- list()
		network_walktrap <- walktrap.community(network_igraph,steps=2)
		modularityMax <- max(network_walktrap$modularity)
		optimalwalktrap <- network_walktrap
		optimalStep <- 2
		
		for(i in c(3:maxStep)){
			network_walktrap <- walktrap.community(network_igraph,steps=i)
			network_modularity <- max(network_walktrap$modularity)
#cat("Modularity:",network_modularity,"\n")
			if(network_modularity>modularityMax){
				optimalwalktrap <- network_walktrap
				optimalStep <- i
				modularityMax <- network_modularity
			}
		}
		
		network_adjMatrix <- calculateRandomWalkerAdjectMatrix(network_igraph,optimalStep)
		network_adjMatrix <- as.dist(network_adjMatrix)
		network_hclust <- transformFromWalktrapToHclust(optimalwalktrap)
		network_order <- seriate(network_adjMatrix,method="OLO",control=list(hclust=network_hclust))
		network_order <- get_order(network_order)
		maxWalktrap <- list(walktrap=optimalwalktrap,step=optimalStep,network=network_igraph,order=network_order,level=level)
		return(maxWalktrap)
	}
	
	identifySig <- function(network_info,method,ModularityThr, ZRandomNum, permuteNum, pThr){
		
	#identify whether the network can be separated again
		
		network_walktrap <- network_info$walktrap
		network_modularity <- max(network_walktrap$modularity)
		network_igraph <- network_info$network
		degree <- igraph:::degree(network_igraph)
		ranmodu <- vector()
		step <- network_info$step
		
		sig <- 0
		
		if(method=="Modularity Cutoff"){
			if(network_modularity>ModularityThr){
				sig <- 1
			}
		}
		
		if(method=="ZScore"){
			for(i in c(1:ZRandomNum)){
				suppressWarnings(rannet <- degree.sequence.game(degree,method="vl"))
				ran_walktrap <- walktrap.community(rannet,steps=step)
				ranModularity <- max(ran_walktrap$modularity)
				ranmodu <- c(ranmodu,ranModularity)
			}
			ranmodu_mean <- mean(ranmodu)
			ranmodu_sd <- sd(ranmodu)
			if(ranmodu_sd == 0 && network_modularity>ranmodu_mean){
				sig <- 1
			}else{
				Z_network_modularity <- (network_modularity-ranmodu_mean)/ranmodu_sd
				p <- pnorm(Z_network_modularity,mean=0,sd=1,lower.tail=F)
				if(p<pThr){
					sig <- 1
				}
			}
		}
		
		if(method=="Permutation"){
			for(i in c(1:permuteNum)){
				suppressWarnings(rannet <- degree.sequence.game(degree,method="vl"))
				ran_walktrap <- walktrap.community(rannet,steps=step)
				ranModularity <- max(ran_walktrap$modularity)
				ranmodu <- c(ranmodu,ranModularity)
			}
			
			p <- length(ranmodu[ranmodu>=network_modularity])/permuteNum
			if(p<pThr){
				sig <- 1
			}
		}
		
		return(sig)
	}
	
	

	identifyHierOr <- function(network_maxComponent,minModule, maxStep, method, ModularityThr, ZRandomNum, permuteNum, pThr){
		
	#identify hierarchical organization of network
		
		allHir <- list()
		ai <- 1
		sigHir <- list()
		si <- 1
		start <- 1
		
		cat("Evaluate Leve 1 network...\n")
		
		levelid <- 2
		
		network_info <- evaluateWalktrapStep(network_maxComponent,maxStep,level=1)
		network_sig <- identifySig(network_info,method,ModularityThr, ZRandomNum, permuteNum, pThr)
		
		allHir[[ai]] <- network_info
		ai <- ai+1
		
		if(network_sig==1){
			sigHir[[si]] <- network_info
			si <- si+1
		}
		
		while(length(sigHir)>0 && start!=si){
			currentNetwork <- sigHir[[start]]
			start <- start+1
			currentNetwork_level <- currentNetwork$level
			currentNetwork_walktrap <- currentNetwork$walktrap
			currentNetwork_igraph <- currentNetwork$network
			currentNetwork_node <- currentNetwork_walktrap$names
			currentNetwork_membership <- currentNetwork_walktrap$membership
			currentNetwork_membership <- data.frame(node=currentNetwork_node,membership=currentNetwork_membership,stringsAsFactors=F)
			currentNetwork_membership_count <- tapply(currentNetwork_membership[,1],currentNetwork_membership[,2],length)
			currentNetwork_membership_count <- currentNetwork_membership_count[currentNetwork_membership_count>=minModule]
			
			if(length(currentNetwork_membership_count)==0){
				next
			}
			
			currentNetwork_membership_group <- as.integer(names(currentNetwork_membership_count))
			
			subLevel <- currentNetwork_level+1
			if(subLevel==levelid){
				cat("Evaluate Level ",levelid," networks...\n",sep="")
				levelid <- levelid+1
			}
			for(i in c(1:length(currentNetwork_membership_group))){
				subnetwork_node <- currentNetwork_membership[currentNetwork_membership[,2]==currentNetwork_membership_group[i],1]
				subnetwork_igraph <- induced.subgraph(currentNetwork_igraph,subnetwork_node)
				subnetwork_info <- evaluateWalktrapStep(subnetwork_igraph,maxStep,subLevel)
				subnetwork_sig <- identifySig(subnetwork_info,method,ModularityThr, ZRandomNum, permuteNum, pThr)
				allHir[[ai]] <- subnetwork_info
				ai <- ai+1
				if(subnetwork_sig==1){
					sigHir[[si]] <- subnetwork_info
					si <- si+1
				}
			}
		}
		
		return(allHir)
	}
	
	
	orderDiffLevel <- function(subnetworkInfo,geneorder){
		
	#reorder all genes in the network according to the optimal position in each level
		
		for(l in c(1:length(subnetworkInfo))){
			allHir <- subnetworkInfo[[l]]
		
			for(i in c(1:length(allHir))){
				ori <- allHir[[i]]$order
				walktrap <- allHir[[i]]$walktrap
				node <- walktrap$names
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
	
	createHMIFile <- function(subnetworkInfo,geneorder){
	#create HMI file
		hmiFile <- data.frame(best="N",level=0,order=1,name="ALL",start=1,end=nrow(geneorder),stringsAsFactors=F)
		hi <- 2
		for(l in c(1:length(subnetworkInfo))){
			allHir <- subnetworkInfo[[l]]
			for(i in c(1:length(allHir))){
				subnetwork <- allHir[[i]]$network
				node <- V(subnetwork)$name
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
	
	
	createGMTFile <- function(hmiFile,geneorder){
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
	
	
	if(missing(inputNetwork)){
		stop("Please input the network!\n")
	}
	
	if(missing(outputFileName)){
		stop("Please input the output file name!\n")
	}else{
		if(substr(outputFileName,nchar(outputFileName),nchar(outputFileName))=="/"){
			stop("Please input the output file name!\n")
		}
	}
	
	#require(igraph) || stop("Package igraph version 0.6 is required!")
	#require(seriation) || stop("Package seriation version 1.0-10 is required!")
	#require(graph) || stop("Package graph version 1.34.0 is required!")
	
	
	
	if(length(which(method %in% c("Modularity Cutoff","ZScore","Permutation")))==0){
		stop("The inputted 'method' is invalid! Please select a method from 'Modularity Cutoff','ZScore' and 'Permutation'!\n")
	}
	
	#load Network
	if(class(inputNetwork)=="character"){
		network <- read.graph(inputNetwork,format="ncol")
	}else{
		if(class(inputNetwork)=="data.frame"){
			if(ncol(inputNetwork)!=2){
				stop("data object should contain two columns!\n");
			}else{
				inputNetwork <- as.matrix(inputNetwork)
				inputNetwork_S <- as.character(inputNetwork)
				inputNetwork_S <- array(inputNetwork_S,dim=dim(inputNetwork))
				network <- graph.edgelist(inputNetwork_S,directed=F)
				rm(inputNetwork,inputNetwork_S)
				gc()
			}
		}else{
			if(class(inputNetwork)=="matrix"){
				if(ncol(inputNetwork)!=2){
					stop("data object should contain two columns!\n");
				}else{
					inputNetwork_S <- as.character(inputNetwork)
					inputNetwork_S <- array(inputNetwork_S,dim=dim(inputNetwork))
					network <- graph.edgelist(inputNetwork_S,directed=F)
					rm(inputNetwork,inputNetwork_S)
					gc()
				}
			}else{
				if(class(inputNetwork)=="graphNEL"){
					network <- igraph.from.graphNEL(inputNetwork)
					rm(inputNetwork)
					gc()
				}else{
					stop("The input network should be from a file or a data object with data.frame or matrix class. Other type of input is invalid!\n")
				}
			}
		}
	}

	proteinInNetwork1 <- V(network)$name
	cat("Network has ",vcount(network)," nodes and ",ecount(network)," edges\n",sep="")
	
	
	network <- simplify(network)
	proteinInNetwork <- V(network)$name
	if(length(proteinInNetwork1) != length(proteinInNetwork)){
		cat("After removing self interactions and loop interactions, network remains ",vcount(network)," nodes and ",ecount(network)," edges\n",sep="")
	}
	
	if(minModule == (-1)){
		if(round(length(proteinInNetwork)*(0.003))>5){
			minModule <- round(length(proteinInNetwork)*(0.003))
		}else{
			minModule <- 5
		}
	}
	
	#if the outputfile is used for NetGestalt (NetGestalt is TRUE), the inputfile should only contain official gene symbols
	
	#if(NetGestalt==TRUE){
		
		#require(org.Hs.eg.db) || stop("Package org.Hs.eg.dbs version 2.8.0 is required!")
		#genesymbolPath <- system.file("extdata","Hsapiens_gene_symbol.txt",package="NetSAM")
		#cat(genesymbolPath,"\n")
		#genesymbol <- read.table(genesymbolPath,header=F,sep="\t",stringsAsFactors=F)
		#cat("all genesymbol:",nrow(genesymbol),"\n")
		#genesymbol <- as.vector(as.matrix(genesymbol))
		#genesymbol <- keys(org.Hs.egSYMBOL2EG)
		
		#overlap_genesymbol_networkP <- intersect(proteinInNetwork, genesymbol)
		#if(length(overlap_genesymbol_networkP) != length(proteinInNetwork)){
        #if(length(overlap_genesymbol_networkP)>(0.8*length(proteinInNetwork))){
        #	network <- induced.subgraph(network,overlap_genesymbol_networkP)
        #	cat("After removing IDs which are not official HUGO Symbols, Network has ",vcount(network)," nodes and ",ecount(network)," edges\n",sep="")
		#	}else{
		#		stop("Over 20% IDs in the inputted network are not official Human HUGO Symbols. Please first transform IDs to official HUGO Symbols and then run NetSAM!\n")
		#	}
		#}
        #}
	
	cat("\nIdentifying the hierarchical modules of the network...\n\n")
	overlap_genesymbol_networkP <- V(network)$name
	subnetworkInfo <- list()
	subnetworkInfo_index <- 1
	network_cluster <- clusters(network)
	network_cluster_size <- network_cluster$csize
	
	if(length(network_cluster_size[network_cluster_size>=minModule])==0){
		stop("The size of all subnetworks in the inputted network are less than ",minModule,". Please adjust the parameter 'minModule'!\n\n")
	}
	
	network_cluster_size <- data.frame(id=c(1:length(network_cluster_size)),cluster_size=network_cluster_size,stringsAsFactors=F)
	network_cluster_size <- network_cluster_size[order(-network_cluster_size[,2]),]
	network_cluster_membership <- network_cluster$membership
	
	network_ordered_node <- vector()
	subnetwork_id <- 1
	for(i in c(1:nrow(network_cluster_size))){
		sub_network_size <- network_cluster_size[i,2]
		if(sub_network_size>=minModule){
			cat("Start to analysis subnetwork ",subnetwork_id,"!\n")
			subnetwork_id <- subnetwork_id+1
			subnetwork_node <- overlap_genesymbol_networkP[which(network_cluster_membership==network_cluster_size[i,1])]
			subnetwork <- induced.subgraph(network,subnetwork_node)
			allHir <- identifyHierOr(subnetwork,minModule, maxStep, method, ModularityThr, ZRandomNum, permuteNum, pThr)
			subnetworkInfo[[subnetworkInfo_index]] <- allHir
			subnetworkInfo_index <- subnetworkInfo_index+1
			network_ordered_node <- c(network_ordered_node,subnetwork_node)
			cat("\n")
		}
		
	}
			
	
	cat("\nReorder the genes in the one dimentional layout...\n")
	
	geneorder <- data.frame(ruler_id=c(1:length(overlap_genesymbol_networkP)),node_type="Gene",node_db="Entrez Gene",node_db_id=overlap_genesymbol_networkP,node_name=overlap_genesymbol_networkP,stringsAsFactors=F)
	rownames(geneorder) <- overlap_genesymbol_networkP
	unann_node <- setdiff(overlap_genesymbol_networkP,network_ordered_node)
	if(length(unann_node) !=0){
		network_ordered_node <- c(network_ordered_node,unann_node)
	}
	geneorder <- geneorder[network_ordered_node,]
	
	
	geneorder <- orderDiffLevel(subnetworkInfo,geneorder)
	
	
	hmiFile <- createHMIFile(subnetworkInfo,geneorder)
	
	#if(NetGestalt==FALSE){
	#	gmtFile <- createGMTFile(hmiFile,geneorder)
	#}
	
	network_edges <- get.edgelist(network)
	
	#if(NetGestalt==TRUE){
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
		
		netgestalt <- list(rulfile=geneorder,hmifile=hmiFile,network=network_edges)

		cat("Processing completed!\n\n") 
		return(netgestalt)
        #}else{
		#outputFile1 <- paste(outputFileName,".gmt",sep="")
		#write.table(gmtFile,file=outputFile1,row.names=F,col.names=F,sep="\t",quote=F)
		#outputFile2 <- paste(outputFileName,".rul",sep="")
		#write.table(geneorder,file=outputFile2,row.names=F,col.names=T,sep="\t",quote=F)
		#outputFile3 <- paste(outputFileName,".net",sep="")
		#write.table(network_edges,file=outputFile3,row.names=F,col.names=T,sep="\t",quote=F)
		#result <- list(gmtfile=gmtFile,rul=geneorder,network=network_edges)
		#cat("Processing completed!\n\n")
		#return(result)
        #}
}

