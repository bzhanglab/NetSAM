mapToSymbol <- function(inputData,organism="hsapiens",inputType="genelist",idType="auto",edgeType="unweighted", collapse_mode="maxSD", is_outputFile=FALSE, outputFileName="", verbose=TRUE){
    
    organisms <- c("hsapiens","mmusculus","rnorvegicus","drerio","celegans","scerevisiae","cfamiliaris","dmelanogaster","athaliana")
    
    if(length(which(organisms==organism))==0){
        stop("Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana!")
    }

    
    if(length(which(inputType %in% c("genelist","network","sbt","sct","matrix")))==0){
        stop("The input 'inputType' is invalide! Please select an option from 'genelist', 'netwrok', 'sbt', 'sct', and 'matrix'!\n")
    }
    
    if(length(which(collapse_mode %in% c("mean","median","maxSD","maxIQR","max","min")))==0){
        stop("The input 'collapse_mode' is invalide! Please select an option from 'mean', 'median', 'maxSD',  'maxIQR', 'max' and 'min' (use mean, median, max standard derivation, max interquartile range, maximum and minimum of duplicate genes for each sample)!\n")
    }
    
    if(inputType=="genelist"){
        
        if(is_outputFile==TRUE && outputFileName==""){
            outputFileName <- "geneList_symbol.txt"
        }
        
        mapresult <- .GeneToSymbol(geneList=inputData,organism=organism,idType=idType,verbose=verbose)
        
        if(is.null(mapresult)){
            if(verbose==TRUE){
                stop("The ids can not be transformed to gene symbols! Please check whether the inut idType is correct!\n",sep="")
            }
            
        }else{
            if(class(mapresult)=="numeric"){
                if(verbose==TRUE){
                    stop(paste("Only ",round(mapresult*100),"% ids can be transformed to gene symbols. Please check whether the inut idType is correct!\n",sep=""))
                }
                mapresult <- NULL
            }else{
                if(is_outputFile==TRUE){
                    write.table(mapresult,file=outputFileName,row.names=F,col.names=F,sep="\t",quote=F)
                }
            }
        }
    }
    
    if(inputType=="network"){
        
        if(is_outputFile==TRUE && outputFileName==""){
            outputFileName <- "network_symbol.net"
        }

        mapresult <- .NetToSymbol(inputNetwork=inputData, organism=organism, idType=idType, edgeType=edgeType, verbose=verbose)
        
        if(is.null(mapresult)){
            if(verbose==TRUE){
                stop("The ids in the input network can not be transformed to gene symbols! Please check whether the inut idType is correct!\n",sep="")
            }

        }else{
            if(class(mapresult)=="numeric"){
                if(verbose==TRUE){
                    stop(paste("Only ",round(mapresult*100),"% ids in the input network can be transformed to gene symbols. Please check whether the inut idType is correct!\n",sep=""))
                }
                mapresult <- NULL
            }else{
                if(is_outputFile==TRUE){
                    net <- mapresult$network
                    write.table(net,file=outputFileName,row.names=F,col.names=F,sep="\t",quote=F)
                }
            }
        }
    }
    
    if(inputType=="matrix"){
        
        if(is_outputFile==TRUE && outputFileName==""){
            outputFileName <- "matrix_symbol.cct"
        }
        
        mapresult <- .MatToSymbol(inputMat=inputData, organism=organism, idType=idType, collapse_mode=collapse_mode, verbose=verbose)
        
        if(is.null(mapresult)){
            if(verbose==TRUE){
                stop("The ids in the input matrix can not be transformed to gene symbols! Please check whether the inut idType is correct!\n")
            }
        }else{
            if(class(mapresult)=="numeric"){
                if(verbose==TRUE){
                    stop(paste("Only ",round(mapresult*100),"% ids in the input matrix can be transformed to gene symbols. Please check whether the inut idType is correct!\n",sep=""))
                }
                mapresult <- NULL
            }else{
                if(is_outputFile==TRUE){
                    matrix <- mapresult$data_symbol
                    matrix <- cbind(rownames(matrix),matrix)
                    colnames(matrix)[1] <- "GeneSymbol"
                    write.table(matrix,file=outputFileName,row.names=F,col.names=T,sep="\t",quote=F)
                }
            }
        }
    }
    
    if(inputType=="sbt"){
        if(is_outputFile==TRUE && outputFileName==""){
            outputFileName <- "sbt_symbol.sbt"
        }

        mapresult <- .sbtToSymbol(inputFile=inputData, organism=organism, idType=idType, verbose=verbose)
        
        if(is.null(mapresult)){
            if(verbose==TRUE){
                stop("The ids in the SBT file can not be transformed to gene symbols! Please check whether the inut idType is correct!\n")
            }
        }else{
            l <- mapresult$oriLength
            mapresult <- mapresult$transformedData
            
            if(l!=nrow(mapresult)){
                if(verbose==TRUE){
                    warning(paste("The ids in ",nrow(mapresult)," tracks among all ",l," tracks in the input SBT file can be transformed to gene symbols!\n",sep=""))
                }
            }
            
            if(is_outputFile==TRUE){
                    write.table(mapresult,file=outputFileName,row.names=F,col.names=F,sep="\t",quote=F)
            }
        }
    }
    
    if(inputType=="sct"){
        
        if(is_outputFile==TRUE && outputFileName==""){
            outputFileName <- "sct_symbol.sct"
        }

        
        mapresult <- .sctToSymbol(inputFile=inputData, organism=organism, idType=idType, verbose=verbose, collapse_mode=collapse_mode)
        
        if(is.null(mapresult)){
            if(verbose==TRUE){
                stop("The ids in the SCT file can not be transformed to gene symbols! Please check whether the inut idType is correct!\n")
            }
        }else{
            if(class(mapresult)=="numeric"){
                if(verbose==TRUE){
                    stop(paste("Only ",round(mapresult*100),"% ids in the SCT file can be transformed to gene symbols. Please check whether the inut idType is correct!\n",sep=""))
                }
                mapresult <- NULL
            }else{
            		sctM <- mapresult$sct_symbol
            		sctM <- cbind(rownames(sctM),sctM)
            		colnames(sctM)[1] <- "GeneSymbol"
                if(is_outputFile==TRUE){
                    write.table(sctM,file=outputFileName,row.names=F,col.names=T,sep="\t",quote=F)
                }
            }
        }
    }
    
    return(mapresult)
}

.GeneToSymbol <-
function(geneList, organism="hsapiens", idType="auto",verbose=TRUE){
	
    
    
    #require(biomaRt) || stop("Package biomaRt version 2.18.0 is required!")
    
    if(!is.null(dim(geneList))){
        stop("The input gene list should be a vector!!\n")
    }
   
   if(organism=="athaliana"){
       mart <- useMart("plants_mart",host="plants.ensembl.org")
       m <- try(mart <- useDataset(paste(organism,"_eg_gene",sep=""),mart),silent=TRUE)
   }else{
       mart <- useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org")
       m <- try(mart <- useDataset(paste(organism,"_gene_ensembl",sep=""),mart),silent=TRUE)
   }
   
    
    if(class(m)=="try-error"){
        stop(paste("The function can not connect to Biomart. Please try again!\n",sep=""))
    }
    
    
    allMartAttr <- listAttributes(mart)
    allFilter <- listFilters(mart)
    
    finalAttr <- intersect(allMartAttr[,1],allFilter[,1])
    
    genesymbol <- ""
    
    if(organism=="hsapiens"){
        noAttr <- c("strand","transcript_status","atlas_diseasestate","transcript_count","phenotype_description","source","go_id","transcript_source","goslim_goa_accession","ens_hs_gene","chromosome_name","status","mim_morbid_accession","atlas_celltype","germ_line_variation_source")
        genesymbol <- "hgnc_symbol"
    }
    
    if(organism=="mmusculus"){
        noAttr <- c("chromosome_name","go_id","goslim_goa_accession","strand","ikmcs_no_products_available_yet","status","germ_line_variation_source","hgnc_id","hgnc_symbol","hgnc_transcript_name","prints","phenotype_description","source")
        genesymbol <- "mgi_symbol"
    }
    
    if(organism=="rnorvegicus"){
        noAttr <- c("chromosome_name","go_id","goslim_goa_accession","strand","status","germ_line_variation_source","mgi_id","mgi_symbol","mgi_transcript_name","prints","source","transcript_source","transcript_count","transcript_biotype","transcript_status")
        genesymbol <- "rgd_symbol"
    }
    
    if(organism=="drerio"){
        noAttr <- c("chromosome_name","go_id","goslim_goa_accession","strand","status","status","germ_line_variation_source","hgnc_id","hgnc_symbol","hgnc_transcript_name","prints","phenotype_description","source","transcript_source","transcript_count","transcript_biotype","transcript_status")
        genesymbol <- "zfin_symbol"
    }
    
    if(organism=="celegans"){
         noAttr <- c("chromosome_name","go_id","prints","source","status","strand","transcript_source","transcript_count","transcript_biotype","transcript_status")
         genesymbol <- "wormbase_locus"
    }
    
    if(organism=="scerevisiae"){
        noAttr <- c("chromosome_name","germ_line_variation_source","go_id","goslim_goa_accession","prints","source","status","strand","transcript_source","transcript_count","transcript_biotype","transcript_status")
        genesymbol <- "external_gene_name"
    }
    
    if(organism=="cfamiliaris"){
        noAttr <- c("chromosome_name","germ_line_variation_source","go_id","goslim_goa_accession","hgnc_id","hgnc_symbol","hgnc_transcript_name","phenotype_description","prints","source","status","strand","transcript_source","transcript_count","transcript_biotype","transcript_status")
        genesymbol <- "external_gene_name"
    }
    
    if(organism=="dmelanogaster"){
        noAttr <- c("chromosome_name","germ_line_variation_source","go_id","goslim_goa_accession","prints","source","status","strand","transcript_source","transcript_count","transcript_biotype","transcript_status")
        genesymbol <- "external_gene_name"
    }
    
    if(organism=="athaliana"){
        noAttr <- c("chromosome_name","prints_pf","strand")
        genesymbol <- "tair_symbol"

    }
    
    finalAttr <- setdiff(finalAttr,noAttr)
    
    finalAttr <- sort(finalAttr)
    
    if(idType!="auto" && length(which(finalAttr==idType))==0){
        cat("The input 'idType' is not supported by Biomart. Please select one id type from the following list:\n")
        cat(finalAttr,sep="\n")
        stop("The users can also set 'idType' as 'auto' and the function will automatically search idtype based on the input data. However, this function may take 20-60 minutes based on users' internet speed\n")
    }
    
    if(idType=="auto"){
        idtype <- .autoSearchIdType(geneList,finalAttr,mart,genesymbol,verbose)
        if(!is.null(idtype)){
            idmap <- getBM(attributes=c(idtype,genesymbol),filters=idtype,values=geneList,mart=mart)
            idmap <- idmap[idmap[,2]!="",]
            
            t <- tapply(idmap[,2],idmap[,1],paste,collapse="|")
            
            idmap <- data.frame(id=names(t),genesymbol=t,stringsAsFactors=F)
            
            return(idmap)
        }else{
            return(NULL)
        }
    }else{
        idmap <- getBM(attributes=c(idType,genesymbol),filters=idType,values=geneList,mart=mart)
        
        if(nrow(idmap)==0){
            return(NULL)
        }else{
            idmap <- idmap[idmap[,2]!="",]
        
            t <- tapply(idmap[,2],idmap[,1],paste,collapse="|")
        
            idmap <- data.frame(id=names(t),genesymbol=t,stringsAsFactors=F)
            
            if(nrow(idmap)/length(geneList)<0.1){
                return(nrow(idmap)/length(geneList))
            }else{
                return(idmap)
            }
        }
    }
}

.autoSearchIdType <- function(geneList,finalAttr,mart,genesymbol,verbose=TRUE){
    if(length(geneList)>40){
        randomNum <- 40
    }else{
        randomNum <- length(geneList)
    }
    
    if(verbose==TRUE){
        cat("\nSearching the array type...\n")
    }
    
    randomId <- sample(geneList,randomNum)
    
    
    finalIdtype <- c()
    tempIdtype <- data.frame(id="",rate=0,stringsAsFactors=F)
    ti <- 1
    
    rang <- ifelse(length(finalAttr)>10,10,length(finalAttr))
    
    startA <- 1
    endA <- rang
    
    startB <- length(finalAttr)-rang+1
    endB <- length(finalAttr)
    
    
    l <- 0
    
    
    u <- 0
    fi <- 0
    while(u<length(finalAttr)){
        if(l==0){
            s <- startA
            e <- endA
        }else{
            s <- startB
            e <- endB
        }
        for(i in c(s:e)){
            idtype <- finalAttr[i]
            idmap <- getBM(attributes=c(idtype,genesymbol),filters=idtype,values=randomId,mart=mart)
            if(length(unique(idmap[,1]))>(0.8*length(randomId))){
                finalIdtype <- idtype
                fi <- 1
                break
            }else{
                if(length(unique(idmap[,1]))>(0.6*length(randomId))){
                    tempIdtype[ti,1] <- idtype
                    tempIdtype[ti,2] <- (length(unique(idmap[,1]))/length(randomId))
                    ti <- ti + 1
                }
            }
        }
        if(fi==1){
            break
        }
        
        if(l==0){
            l <- 1
            u <- u+endA-startA+1
            startA <- startA + rang
            endA <- startA + rang - 1
            if(endA>=startB){
                endA <- startB-1
            }
        }else{
            l <- 0
            u <- u+endB-startB+1
            endB <- startB-1
            startB <- endB - rang + 1
            if(startB<endA){
                startB <- endA+1
            }
        }
        
    }
    
    is.map <- TRUE
    if(length(finalIdtype)>0){
        idtype <- finalIdtype
        is.map <- TRUE
        
    }else{
        if(tempIdtype[1,1]!=""){
            tempIdtype <- tempIdtype[order(-tempIdtype[,2]),]
            idtype <- tempIdtype[1,1]
            is.map <- TRUE
        }else{
            is.map <- FALSE
        }
    }
    
    if(is.map==TRUE){
        
        if(verbose==TRUE){
            cat("The input Ids can map to ",idtype, "!\n")
        }
        
        return(idtype)
        
    }else{
        if(verbose==TRUE){
            cat("The function can not find the matched id type from biomart!\n")
        }
        return(NULL)
    }
}


.NetToSymbol <-
function(inputNetwork, organism="hsapiens", idType="auto", edgeType="unweighted", verbose=TRUE){
    
    
    #require(tools) || stop("Package tools version 3.0.0 is required!")
    
    #load Network
    
    if(verbose==TRUE){
        cat("Transforming the ids in the input network to gene symbols...\n")
    }
    
    
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
    
    netNode <- get.vertex.attribute(network,"name")
    idmap <- .GeneToSymbol(netNode,organism,idType,verbose)
    
    if(is.null(idmap)){
        return(NULL)
    }else{
        if(class(idmap)=="numeric"){
            return(idmap)
        }else{
            network <- .extractSubNetwork(network,idmap[,1],edgeType)
            rownames(idmap) <- idmap[,1]
            idmap <- idmap[get.vertex.attribute(network,"name"),]
            V(network)$name <- idmap[,2]
            network <- simplify(network)
            # if(verbose==TRUE){
            #    cat("After mapping to the gene symbol, the network contains ",vcount(network)," nodes and ",ecount(network)," edges!\n",sep="")
            #}
        
            network_edgelist <- get.edgelist(network)
            network_edgelist <- as.data.frame(network_edgelist)
            if(edgeType=="weighted"){
                wei <- get.edge.attribute(network,name="weight")
                network_edgelist$weight <- wei
            }
        
            re <- list(network=network_edgelist,idmap=idmap)
            return(re)
        }
    }
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


.MatToSymbol <-
function(inputMat, organism="hsapiens", idType="auto", collapse_mode="maxSD", verbose=TRUE){
    
    
    #require(biomaRt) || stop("Package biomaRt version 2.18.0 is required!")
    #require(tools) || stop("Package tools version 3.0.0 is required!")
    #require(multicore) || stop("Package multicore version 0.1-7 is required!")
    #require(doSNOW) || stop("Package doSNOW version 1.0.6 is required!")
    #require(foreach) || stop("Package foreach version 1.4.0 is required!")
    
    #mart <- useMart("ensembl")
    
    #m <- try(mart <- useDataset(paste(organism,"_gene_ensembl",sep=""),mart),silent=TRUE)
    #if(class(m)=="try-error"){
    #    stop(paste("The organism ",organism," your input is not valid! Correct organism names can be obtained with the listDatasets function.\n",sep=""))
    #}

    
    re <- testFileFormat(inputMat=inputMat, collapse_mode=collapse_mode)
    inputMat <- re$inputMat
    inputId <- rownames(inputMat)
    
    idmap <- .GeneToSymbol(inputId,organism,idType,verbose)
    if(!is.null(idmap)){
        
        if(class(idmap)=="numeric"){
            return(idmap)
        }else{
            inputMat <- inputMat[idmap[,1],]
            id <- as.vector(idmap[,2])
            inputMat <- mergeDuplicate(id,inputMat,collapse_mode)
            re <- list(data_symbol=inputMat,idmap=idmap)
            return(re)
        }
    }else{
        return(NULL)
    }
}


.sbtToSymbol <- function(inputFile, organism="hsapiens", idType="auto", verbose=TRUE){
    if(class(inputFile)=="character"){
        if(file_ext(inputFile)!="sbt"){
            stop("The extension of the input file should be 'sbt'. The detail of the 'sbt' file format can be found in the NetGestalt (www.netgestalt.org)!\n")
        }else{
            inputMat <- read.csv(inputFile,header=FALSE,sep="\t",stringsAsFactors=FALSE)
        }
    }else{
        stop("The 'inputFile' should be the directory of a SBT file\n")
    }
    
    .calculateLength <- function(vector){
        vector <- vector[vector!=""]
        return(length(vector))
    }
    
    sbtLength <- apply(inputMat,1,.calculateLength)
    
    if(length(which(sbtLength<3))>0){
        stop("The SBT file should contain at least three columns. The following rows of the input file only contain one or two columns: ",which(sbtLength<3),"\n",sep="")
    }
    
    transformedData <- data.frame(name="",des="",genelist="",stringsAsFactors=F)
    ti <- 1
    for(i in c(1:nrow(inputMat))){
        geneL <- inputMat[i,3:ncol(inputMat)]
        geneL <- geneL[geneL!=""]
        idmap <- .GeneToSymbol(geneL,organism,idType,verbose)
        mapNum <- c()
        if(!is.null(idmap) && class(idmap)!="numeric"){
            mappedL <- unique(idmap[,2])
            transformedData[ti,1] <- inputMat[i,1]
            transformedData[ti,2] <- inputMat[i,2]
            mappedL <- paste(mappedL,collapse="\t")
            transformedData[ti,3] <- mappedL
            ti <- ti+1
        }
    }
    
    if(transformedData[1,1]!=""){
        re <- list(transformedData=transformedData,oriLength=nrow(inputMat))
        return(re)
    }else{
        return(NULL)
    }
}


.sctToSymbol <- function(inputFile, organism="hsapiens", idType="auto", verbose=TRUE, collapse_mode="min"){
    if(class(inputFile)=="character"){
        if(file_ext(inputFile)!="sct"){
            stop("The extension of the input file should be 'sct'. The detail of the 'sct' file format can be found in the NetGestalt (www.netgestalt.org)!\n")
        }else{
            inputMat <- read.csv(inputFile,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }
    }else{
        stop("The 'inputFile' should be the directory of a SCT file\n")
    }
    
    if(ncol(inputMat)<2){
        stop("The SCT file should contain at least two columns (gene id and statistic score). \n")
    }
    
    id <- inputMat[,1]
    inputMat <- as.matrix(inputMat[,2:ncol(inputMat)])
    
    if(length(unique(id))<length(id)){
        cat("The SCT file contains duplidate Ids. The function will use 'median' to collapse duplicate Ids for each statistic.\n")
        inputMat <- mergeDuplicate(id,inputMat,collapse_mode)
    }else{
        rownames(inputMat) <- id
    }
    
    
    idmap <- .GeneToSymbol(id,organism,idType,verbose)
    
    if(!is.null(idmap)){
        
        if(class(idmap)=="numeric"){
            return(idmap)
        }else{
            inputMat <- as.matrix(inputMat[idmap[,1],])
            id <- as.vector(idmap[,2])
            inputMat <- mergeDuplicate(id,inputMat,collapse_mode)
            re <- list(sct_symbol=inputMat,idmap=idmap)
            return(re)
        }
    }else{
        return(NULL)
    }
}


