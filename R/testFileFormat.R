testFileFormat <- function(inputMat=NULL,sampleAnn=NULL,collapse_mode="maxSD"){
    
    re <- list(inputMat=NULL,sampleAnn=NULL)
    
    if(length(which(collapse_mode %in% c("mean","median","maxSD","maxIQR","max","min")))==0){
        stop("The input 'collapse_mode' is not valide! Please select an option from 'mean', 'median', 'maxSD',  'maxIQR', 'max' and 'min' (use mean, median, max standard derivation, max interquartile range, maximum and minimum of duplicate genes for each sample)!\n")
    }
    
    if(!is.null(inputMat) || !is.null(sampleAnn)){
        if(!is.null(inputMat)){
            inputMat <- .testCCTFormat(inputMat=inputMat,collapse_mode=collapse_mode)
            re$inputMat <- inputMat
        }
        if(!is.null(sampleAnn)){
            sampleAnn <- .testTSIFormat(inputMat=inputMat,sampleAnn=sampleAnn)
            re$sampleAnn <- sampleAnn
        }
        return(re)
    }else{
        stop("Please input at least data matrix or sample annotation data!")
    }
}


.testCCTFormat <- function(inputMat, collapse_mode="maxSD"){
    
    if(class(inputMat)=="character"){
        if(file_ext(inputMat)!="cct" && file_ext(inputMat)!="cbt"){
            stop("The extension of the input file should be 'cct' or 'cbt'. The detail of the 'cct' or 'cbt' file format can be found in the NetGestalt (www.netgestalt.org)!\n")
        }else{
            inputMat <- read.table(inputMat,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
            geneid <- inputMat[,1]
            
            if(length(geneid)!=length(unique(geneid))){
                cat("The input data contain duplicate Id. The function will use ",collapse_mode," to collapse duplicate Id in each sample!\n",sep="")
                inputMat <- mergeDuplicate(geneid,inputMat[,2:ncol(inputMat)] ,collapse_mode)
            }else{
                inputMat <- inputMat[,c(2:ncol(inputMat))]
                rownames(inputMat) <- geneid
            }
        }
    }else{
        if(class(inputMat) != "matrix" && class(inputMat) != "data.frame"){
            stop("The type of input data should be a matrix or data.frame object. Other types of data are not supported by this package.!\n")
        }else{
            x <- apply(inputMat,2,function(e) return(class(e)=="numeric" || class(e)=="integer"))
            y <- all(x==TRUE)
            if(y==FALSE){
                stop("The input matrix or data.frame object should only contain numeric or integer values.\n")
            }
        }
    }
    
    if(ncol(inputMat)<6){
        stop("The data should contain at least six samples!\n")
    }
    
    if(length(which(inputMat %in% c("Inf","-Inf")))>0){
        stop("The input data contain Inf which may be gernated by some wrong operation, such as log(0) or 1/0. Please re-process the data and remove the Inf\n")
    }
    
    return(inputMat)
}

.testTSIFormat <- function(inputMat,sampleAnn){
    
    if(class(sampleAnn)=="character"){
        if(file_ext(sampleAnn)!="tsi"){
            stop("The extension of the input annotation file should be 'tsi'. The detail of the 'tsi' file format can be found in the NetGestalt (www.netgestalt.org)!\n")
        }else{
            sampleAnn <- read.table(sampleAnn,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
        }
    }else{
        if(class(sampleAnn) != "data.frame"){
            stop("The type of input annotation data should be a data.frame object. Other types of data are not supported by this package.!\n")
        }
    }
    
    data_type <- sampleAnn[1,2:ncol(sampleAnn)]
    if(length(setdiff(unique(data_type),c("BIN","CAT","CON","SUR")))>0){
        stop("The second row of the tsi file is the data type of the annotations! The current version supports four types of features: BIN, CAT, CON, SUR. The input annotation file also contains other types:",setdiff(unique(data_type),c("BIN","CAT","CON","SUR")),".\n")
    }
    
    if(sampleAnn[2,1]=="category"){
        data <- sampleAnn[3:nrow(sampleAnn),2:ncol(sampleAnn)]
        annSampleName <- sampleAnn[3:nrow(sampleAnn),1]
        startR <- 3
    }else{
        data <- sampleAnn[2:nrow(sampleAnn),2:ncol(sampleAnn)]
        annSampleName <- sampleAnn[2:nrow(sampleAnn),1]
        startR <- 2
    }
    
    sampleName <- colnames(inputMat)
    sampleName <- data.frame(id=c(1:length(sampleName)),name=sampleName,stringsAsFactors=FALSE)
    sampleName <- sampleName[order(sampleName[,2]),]
    annSampleName <- data.frame(id=c(1:length(annSampleName)),name=annSampleName,stringsAsFactors=FALSE)
    annSampleName <- annSampleName[order(annSampleName[,2]),]
    
    if(sum(sampleName[,2]==annSampleName[,2])!=nrow(sampleName)){
        stop("The sample names in the matrix data and annotation data should be exactly same!")
    }
    
    
    for(i in c(1:length(data_type))){
        dt <- data_type[i]
        d <- data[,i]
        if(dt=="BIN"){
            d <- d[!is.na(d)]
            if(length(unique(d))!=2){
                stop(paste("column ",i+1,": Binary type must have 2 distinct values, current column contains ",length(unique(d))," values!\n",sep=""))
            }
        }
        if(dt=="CAT"){
            d <- d[!is.na(d)]
            if(length(d)<3){
                stop(paste("column ",i+1,": Category type must have at least 3 distinct values, current column contains ",length(d)," distinct values!\n",sep=""))
            }
        }
        if(dt=="CON"){
            d[is.na(d)] <- 0
            d <- as.numeric(d)
            if(length(which(is.na(d)))>0){
                stop(paste("column ",i+1,": Contenous type must only contain numbers. Row ",which(is.na(d))+startR," contain characters!\n",sep=""))
            }
        }
        if(dt=="SUR"){
            if(length(which(is.na(d)))>0){
                stop(paste("column ",i+1,": the format of 'NA' for survival information is 'NA,NA' instead of just 'NA'!. Row ",which(is.na(d))+startR," contain 'NA'.\n",sep=""))
            }
            d <- strsplit(d,",",fixed=TRUE)
            l <- lapply(d,function(e){return(length(e))})
            l <- unlist(l)
            if(length(which(l==1))>0){
                stop(paste("column ",i+1,": the format of the survival information is 'time,event'. Row ",which(l==1)+startR," only contain time or event.\n",sep=""))
            }
            d <- do.call(rbind,d)
            d1 <- d[,1]
            d2 <- d[,2]
            
            d1_in <- which(d1=="NA")
            d2_in <- which(d2=="NA")
            if(suppressWarnings(sum(d1_in==d2_in))!=length(d1_in)){
                stop(paste("column ",i+1,": if the survial time is NA, the corresponding event should be NA. Row ",c(setdiff(d1_in,d2_in),setdiff(d2_in,d1_in))," only contain one NA for time or event!\n",sep=""))
            }
            
            d1[d1=="NA"] <- 10
            d1 <- as.numeric(d1)
            if(length(which(is.na(d1)))>0){
                stop(paste("column ",i+1,": the survival time must only be numeric. Row ",which(is.na(d1))+startR," contain characters!\n",sep=""))
            }
            d2 <- d[,2]
            d2[d2=="NA"] <- 0
            d2 <- as.numeric(d2)
            if(length(which(is.na(d2)))>0){
                stop(paste("column ",i+1,": the survival event must only be 0 or 1. Row ",which(is.na(d2))+startR," contain characters!\n",sep=""))
            }
            d2 <- unique(d2)
            if(length(d2)>2){
                stop(paste("column ",i+1,": the survival event must only be 0 or 1. The current data contain ",length(d2)," values!\n",sep=""))
            }
        }
    }
    
    x <- cbind(sampleName,annSampleName)
    x <- x[order(x[,1]),]
    if(sampleAnn[2,1]=="category"){
        y <- sampleAnn[3:nrow(sampleAnn),]
        y <- y[x[,3],]
        sampleAnn <- rbind(sampleAnn[1:2,],y)
    }else{
        y <- sampleAnn[2:nrow(sampleAnn),]
        y <- y[x[,3],]
        sampleAnn <- rbind(sampleAnn[1,],y)
    }
    colnames(sampleAnn)[1] <- "Barcode"
    
    return(sampleAnn)
}