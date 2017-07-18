MatSAM <- function(inputMat, sampleAnn=NULL, outputFileName, outputFormat="msm", organism="hsapiens", map_to_symbol=FALSE, idType="auto", collapse_mode="maxSD", naPer=0.7, meanPer=0.8, varPer=0.8, corrType="spearman", matNetMethod="rank", valueThr=0.5, rankBest=0.003, networkType="signed", netFDRMethod="BH", netFDRThr=0.05, minModule=0.003, stepIte=FALSE, maxStep=4, moduleSigMethod="cutoff", modularityThr=0.2, ZRanNum=10, PerRanNum=100, ranSig=0.05, idNumThr=(-1), nThreads=3){
    
    #require(tools) || stop("Package tools version 3.0.0 is required!")
    
    organisms <- c("hsapiens","mmusculus","rnorvegicus","drerio","celegans","scerevisiae","cfamiliaris","dmelanogaster", "athaliana")
    names(organisms) <- c("Hs","Mm","Rn","Dr","Ce","Sc","Cf","Dm","At")
    
    if(outputFormat != "none" && missing(outputFileName)){
        stop("Please input the output file name!\n")
    }else{
        if(substr(outputFileName,nchar(outputFileName),nchar(outputFileName))=="/"){
            stop("Please input the output file name!\n")
        }
    }
    
    if(length(which(outputFormat %in% c("msm","gmt","multiple","none")))==0){
        stop("The input 'outputFormat' is invalid! Please select an output format from 'nsm', 'gmt', 'multiple' and 'none'!\n")
    }

    
    if(length(which(organisms==organism))==0){
        stop("Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana!")
    }
    
    if(length(which(collapse_mode %in% c("mean","median","maxSD","maxIQR","max","min")))==0){
        stop("The input 'collapse_mode' is not valide! Please select an option from 'mean', 'median', 'maxSD',  'maxIQR', 'max' and 'min' (use mean, median, max standard derivation, max interquartile range, maximum and minimum of duplicate genes for each sample)!\n")
    }


    if(meanPer>1 || meanPer<0 || class(meanPer)!="numeric"){
        stop("The input 'meanPer' is invalid! 'meanPer' should be a positive number and less than 1!\n")
    }
    
    if(varPer>1 || varPer<0 || class(varPer)!="numeric"){
        stop("The input 'varPer' is invalid! 'varPer' should be a positive number and less than 1!\n")
    }
    
    
    if(length(which(corrType %in% c("pearson","spearman")))==0){
		stop("The input 'corrType' is invalid! Please select a method from 'pearson' and 'spearman'!\n")
	}
    
    if(length(which(matNetMethod %in% c("value","rank")))==0){
		stop("The input 'matNetMethod' is invalid! Please select a method from 'value' and 'rank'!\n")
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
		stop("The input 'netFDRMethod' is invalid! Please select an option from 'holm','hochberg', 'hommel', 'bonferroni', 'BH', 'BY' and 'none'!\n")
	}
    
    if(netFDRThr>1 || class(netFDRThr)!="numeric"){
        stop("The input 'netFDRThr' is invalid! 'netFDRThr' should be a number and less than 1!\n")
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
    
    if(idNumThr<(-1) || class(idNumThr)!="numeric"){
        stop("The input 'idNumThr' is invalid! 'idNumThr' should be -1 or a positive number!\n")
    }
    
    if(nThreads<1 || class(nThreads)!="numeric"){
        stop("The input 'nThreads' is invalid! 'nThreads' should be a number and greater than 1!\n")
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
    ################################Identify whether the input matrix and annotation file have the correct format##############################
    
    #cat("Test file\n")
    
    re <- testFileFormat(inputMat=inputMat,sampleAnn=sampleAnn,collapse_mode=collapse_mode)
    inputMat <- re$inputMat
    
    #cat("Dimension1:",dim(inputMat),"\n")
    
    is_sample <- 0
    if(!is.null(sampleAnn)){
        sampleAnn <- re$sampleAnn
        is_sample <- 1
    }
    
    geneid <- rownames(inputMat)
    sampleName <- colnames(inputMat)
    
    
    ###############Identify whether the IDs in the expression data are gene symbols###############
    
    #cat("Gene Symbol \n")
    
    is.genesymbol <- FALSE
    if(length(intersect(geneid,genesymbol))>0.6*length(geneid)){
        is.genesymbol <- TRUE
    }else{
        is.genesymbol <- FALSE
    }

        
    idmap <- NULL
    mapping.status <- 1
    if(map_to_symbol==TRUE){
        
        re <- mapToSymbol(inputData=inputMat, organism=organism, inputType="matrix", idType=idType, collapse_mode=collapse_mode, verbose=TRUE)
        
        if(!is.null(re)){
            inputMat <- re$data_symbol
            mapping.status <- 0
            idmap <- re$idmap
        }
    }else{
        if(is.genesymbol==TRUE){
            mapping.status <- 0
        }else{
            if(idType==""){      #NONE OF THE ABOVE
                mapping.status <- 1
            }else{
                re <- mapToSymbol(inputData=inputMat, organism=organism, inputType="matrix", idType=idType, collapse_mode=collapse_mode, verbose=FALSE)
             
                if(is.null(re)){
                    stop("The input idType can not match the id type of the input data! Please check whether the inut idType is correct!\n")
                }else{
                    idmap <- re$idmap
                    mapping.status <- 2
                }
            }
        }
    }
    
    ###########################Construct co-expression network##################################
    
    #cat("MatNet\n")
    
    inputMat <- as.matrix(inputMat)
    
    coExp <- MatNet(inputMat=inputMat,collapse_mode=collapse_mode,naPer=naPer,meanPer=meanPer,varPer=varPer,corrType=corrType,matNetMethod=matNetMethod,valueThr=valueThr,rankBest=rankBest,networkType=networkType,netFDRMethod=netFDRMethod,netFDRThr=netFDRThr,idNumThr=idNumThr,nThreads=nThreads)
    filtered_gene <- union(coExp[,1],coExp[,2])
    
    filteredExp <- inputMat[filtered_gene,]
    
    #############################Identify Modules###########################
    
    #cat("NetSAM\n")
    
    netgestalt <- NetSAM(inputNetwork=coExp, outputFileName=outputFileName, outputFormat="none", edgeType="unweighted", map_to_genesymbol=FALSE, organism=organism, idType=idType, minModule=minModule, stepIte=stepIte, maxStep=maxStep, moduleSigMethod=moduleSigMethod, modularityThr=modularityThr, ZRanNum=ZRanNum, PerRanNum=PerRanNum, ranSig=ranSig, edgeThr=(-1), nodeThr=(-1), nThreads=nThreads)
    
    #############################Calcualte relationship between module and sample annotation#########
    hmi <- netgestalt$hmifile
    
    hmi[,7] <- ""
    colnames(hmi)[7] <- "AssociatedFeatures"
    
    if(is_sample==1){
        cat("Calculate the associations between modules and annotations...\n")
        moduleAsso <- featureAssociation(filteredExp,sampleAnn,netgestalt,paste(outputFileName,"_featureAsso",sep=""))
        attrName <- c("FeatureName","Type","StatisticMethod","GroupComparison","Statistics","Direction","PValue","FDR")
    
        if(!is.null(moduleAsso)){
    
            for(i in c(2:nrow(hmi))){
                mN <- hmi[i,4]
                fe <- moduleAsso[moduleAsso[,1]==mN,2:ncol(moduleAsso)]
                fe <- apply(fe,1,.paste_dataframe,attrName,collapse="|")
                fe <- paste(fe,collapse=";")
                hmi[i,7] <- fe
            }
        }
    }
    
    netgestalt$hmifile <- hmi
    
    ########################Identify the most significant GO terms############
    
    if(mapping.status!=1){
        cat("Identify the associated GO Terms for each module...\n")
        goAsso <- GOAssociation(netgestalt,paste(outputFileName,"_GOAsso",sep=""),organism=organism,outputType="top",topNum=1)
        hmi <- netgestalt$hmifile
        
        hmi[,8] <- ""
        colnames(hmi)[8] <- "Related GO Terms"
        
        attrName <- c("Ontology","GOID","GOName","PValue","FDR")
        
        if(!is.null(goAsso)){
            
            for(i in c(2:nrow(hmi))){
                mN <- hmi[i,4]
                fe <- goAsso[goAsso[,1]==mN,2:ncol(goAsso)]
                fe <- apply(fe,1,.paste_dataframe,attrName,collapse="|")
                fe <- paste(fe,collapse=";")
                hmi[i,8] <- fe
            }
        }
        netgestalt$hmifile <- hmi

    }else{
        cat("Because the Ids in the matrix can not be transformed to gene symbol, GO analysis can not be performed!\n")
    }
    
    ##############################Output###################################
    
    d <- dist(t(filteredExp))
    h <- hclust(d)
    filteredExp <- filteredExp[,h$order]
    
    
    exp_output <- cbind(rownames(filteredExp),filteredExp)
    colnames(exp_output)[1] <- "GeneSymbol"
    
    netgestalt$expression <- filteredExp
    
    if(is_sample==1){
        rownames(sampleAnn) <- sampleAnn[,1]
        
        
        if(sampleAnn[2,1]=="category"){
            reo <- c("data_type","category",colnames(filteredExp))
            sampleAnn <- sampleAnn[reo,]
            hierarchicalCluster <- c("CON","cluster",1:ncol(filteredExp))
            sampleAnn <- cbind(sampleAnn,hierarchicalCluster)
        }else{
            reo <- c("data_type",colnames(filteredExp))
            sampleAnn <- sampleAnn[reo,]
            hierarchicalCluster <- c("CON",1:ncol(filteredExp))
            sampleAnn <- cbind(sampleAnn,hierarchicalCluster)
        }
    }
    
    
    if(outputFormat=="msm"){
        outputFile <- paste(outputFileName,".msm",sep="")
        note <- "##  Ruler file  ##"
        write.table(note,file=outputFile,row.names=F,col.names=F,quote=F)
        suppressWarnings(write.table(netgestalt$rulfile,file=outputFile,row.names=F,col.names=T,append=T,quote=F,sep="\t"))
        note <- "##  HMI file ##"
        write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
        suppressWarnings(write.table(netgestalt$hmifile,file=outputFile,row.names=F,col.names=T,quote=F,append=T,sep="\t"))
        note <- "##  Network file  ##"
        write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
        write.table(netgestalt$network,file=outputFile,row.names=F,col.names=F,quote=F,append=T,sep="\t")
    
        note <- "##  Expression data  ##"
        write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
    
                suppressWarnings(write.table(exp_output,file=outputFile,row.names=F,col.names=T,quote=F,append=T,sep="\t"))
    
    
        if(is_sample==1){
            
            note <- "##  Sample Annotation  ##"
            write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
            suppressWarnings(write.table(sampleAnn,file=outputFile,row.names=F,col.names=T,quote=F,append=T,sep="\t"))
            netgestalt$sampleAnn <- sampleAnn
        }
    
        note <- "##  Mapping Status  ##"
        write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
        suppressWarnings(write.table(paste("mapping.status=",mapping.status,sep=""),file=outputFile,row.names=F,col.names=F,quote=F,append=T,sep="\t"))
        netgestalt$mapping.status <- mapping.status
    
        if(mapping.status==2){
            note <- "##  Gene Annotation  ##"
            write.table(note,file=outputFile,row.names=F,col.names=F,quote=F,append=T)
            suppressWarnings(write.table(idmap,file=outputFile,row.names=F,col.names=T,quote=F,append=T,sep="\t"))
            netgestalt$geneAnn <- idmap
        }
    }
    
    if(outputFormat=="multiple"){
        outputFile <- paste(outputFileName,"_rulerFile.rul",sep="")
        write.table(netgestalt$rulfile,file=outputFile,row.names=F,col.names=T,quote=F,sep="\t")
        outputFile <- paste(outputFileName,"_hmiFile.hmi",sep="")
        write.table(netgestalt$hmifile,file=outputFile,row.names=F,col.names=T,quote=F,sep="\t")
        outputFile <- paste(outputFileName,"_corrNet.net",sep="")
        write.table(netgestalt$network,file=outputFile,row.names=F,col.names=F,quote=F,sep="\t")
        outputFile <- paste(outputFileName,"_filteredMatrix.cct",sep="")
        write.table(exp_output,file=outputFile,row.names=F,col.names=T,quote=F,sep="\t")
        if(is_sample==1){
            
            outputFile <- paste(outputFileName,"_standardizedSampleAnn.tsi",sep="")
            write.table(sampleAnn,file=outputFile,row.names=F,col.names=T,quote=F,sep="\t")
            netgestalt$sampleAnn <- sampleAnn
        }
    }
    
    if(outputFormat=="gmt"){
        gmtFile <- .createGMTFile(netgestalt$hmifile,netgestalt$rulfile)
        outputFile1 <- paste(outputFileName,".gmt",sep="")
        write.table(gmtFile,file=outputFile1,row.names=F,col.names=F,sep="\t",quote=F)
        netgestalt$gmt <- gmtFile
    }
    
    return(netgestalt)
    
}


.paste_dataframe <- function(dataframe,attrName,collapse){
    vector <- as.vector(as.matrix(dataframe))
    vector <- paste(attrName,vector,sep=":")
    vector <- paste(vector,collapse=collapse)
    return(vector)
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

