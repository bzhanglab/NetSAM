NetAnalyzer <-
function(inputNetwork, outputFileName, edgeType="unweighted"){
   
  #require(igraph) || stop("Package igraph version 0.6 is required!")
    
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
	
	if(missing(outputFileName)){
		stop("Please input the output file name!\n")
	}else{
		if(substr(outputFileName,nchar(outputFileName),nchar(outputFileName))=="/"){
			stop("Please input the output file name!\n")
		}
	}
	
	nodeAttr <- data.frame(node=vertex_attr(network,name="name"),degree=0,clusterCoefficient=0,betweennessCentrality=0,closenessCentrality=0,stringsAsFactors=F)
	
	if(edgeType=="unweighted"){
	
		nodeAttr[,2] <- degree(network)
		nodeAttr[,3] <- transitivity(network,type="local")
		nodeAttr[,4] <- betweenness(network,directed=FALSE)
		nodeAttr[,5] <- closeness(network,mode="all")
	}else{
		nodeAttr[,2] <- strength(network)
		nodeAttr[,3] <- transitivity(network,type="weighted")
		nodeAttr[,4] <- betweenness(network,directed=FALSE)
		nodeAttr[,5] <- closeness(network,mode="all")
	}
	
	shortestP <- distances(network)
	shO <- paste(outputFileName,"_shortestPathDistance.txt",sep="")
	write.table(shortestP,file=shO,row.names=T,col.names=T,sep="\t",quote=F)
	
	shortestP_V <- shortestP[lower.tri(shortestP)]
	
	nodeO <- paste(outputFileName,".txt",sep="")
	write.table(nodeAttr,file=nodeO,row.names=F,col.names=F,sep="\t",quote=F)
	
	degreeP <- paste(outputFileName,"_degreeDistribution.pdf",sep="")
	if(edgeType=="unweighted"){
		.fit_power_law(network,edgeType,degreeP)
	}else{
		pdf(degreeP,height=6,width=6)
		hist(nodeAttr[,2],xlab="Degree",col="black",main="Degree Distribution")
		dev.off()
	}

	
	shortP <- paste(outputFileName,"_shortestPathDistanceFrequency.pdf",sep="")
	pdf(shortP,height=6,width=6)
	hist(shortestP_V,xlab="Shortest Path Distance",col="black",main="Shortest Path Distance Distribution")
	dev.off()
	
	x <- nodeAttr[nodeAttr[,2]>1,]
	ccP <- paste(outputFileName,"_clusteringCoefficient.pdf",sep="")
	pdf(ccP,height=6,width=6)
	a <- tapply(x[,3],x[,2],mean,na.rm=TRUE)
	plot(as.numeric(names(a)),a,xlab="Degree",ylab="Ave clustering coefficient",main="Clustering Coefficient Distribution",pch=19)
	dev.off()
	

	ccP <- paste(outputFileName,"_betweeness.pdf",sep="")
	pdf(ccP,height=6,width=6)
	a <- tapply(nodeAttr[,4],nodeAttr[,2],mean,na.rm=TRUE)
	plot(as.numeric(names(a)),a,xlab="Degree",ylab="Ave betweeness",main="Betweeness Distribution",pch=19)
	dev.off()
	
	ccP <- paste(outputFileName,"_closenessCentrality.pdf",sep="")
	pdf(ccP,height=6,width=6)
	a <- tapply(nodeAttr[,5],nodeAttr[,2],mean,na.rm=TRUE)
	plot(as.numeric(names(a)),a,xlab="Degree",ylab="Ave closeness centrality",main="Closeness Centrality Distribution",pch=19)
	dev.off()
	
}


.fit_power_law = function(graph,edgeType,outputFile) {
    # calculate degree

    d = degree(graph)
   	dd = degree.distribution(graph, cumulative = FALSE)

    degree = 1:max(d)
    probability = dd[-1]
    # delete blank values
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
    print(paste("Alpha =", round(alpha, 3)))
    print(paste("R square =", round(R.square, 3)))
    # plot
    pdf(outputFile,height=6,width=6)
    plot(probability ~ degree, log = "xy", xlab = "Degree (log)", ylab = "Probability (log)", 
        col = 1, main = "Degree Distribution")
    curve(power.law.fit, col = "red", add = T, n = length(d))
    dev.off()
}
