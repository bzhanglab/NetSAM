featureAssociation <-
function(inputMat,sampleAnn,NetSAMOutput,outputHtmlFile,CONMethod="spearman",CATMethod="kruskal",BINMethod="ranktest",fdrmethod="BH",pth=0.05,collapse_mode="maxSD"){
    
    #require(WGCNA)
    #require(survival)
    #require(R2HTML)
    
   #if(missing(inputMat) || missing(sampleAnn) || sampleAnn==""){
   if(missing(inputMat) || missing(sampleAnn)){
        stop("Please input the data matrix and sample annotation file!\n")
    }

    re <- testFileFormat(inputMat=inputMat, sampleAnn=sampleAnn, collapse_mode=collapse_mode)
    inputMat <- re$inputMat
    sampleAnn <- re$sampleAnn
    
    
    if(length(which(CONMethod %in% c("pearson","spearman")))==0){
        stop("The input 'CONMethod' is invalid! Please select a method from 'pearson' and 'spearman'!\n")
    }
    
    if(length(which(CATMethod %in% c("anova","kruskal")))==0){
        stop("The input 'CATMethod' is invalid! Please select a method from 'anova' and 'kruskal'!\n")
    }
    
    if(length(which(BINMethod %in% c("ttest","ranktest")))==0){
        stop("The input 'BINMethod' is invalid! Please select a method from 'ttest' and 'ranktest'!\n")
    }
    
    if(is.null(NetSAMOutput$rulfile) || is.null(NetSAMOutput$hmifile)){
        stop("NetSAMOutput should be the output from NetSAM function, which contain rulfile, hmifile and network!\n")
    }
    
    if(length(which(collapse_mode %in% c("mean","median","maxSD","maxIQR","max","min")))==0){
        stop("The input 'collapse_mode' is not valide! Please select an option from 'mean', 'median', 'maxSD',  'maxIQR', 'max' and 'min' (use mean, median, max standard derivation, max interquartile range, maximum and minimum of duplicate genes for each sample)!\n")
    }

    
    if(sampleAnn[2,1]=="category"){
        sampleAnn <- sampleAnn[c(1,3:nrow(sampleAnn)),]
    }
    
    

    rule <- NetSAMOutput$rulfile
    hmi <- NetSAMOutput$hmifile
    
    
    if(length(setdiff(rule[,4],rownames(inputMat)))>0){
        stop("All ids in the 'NetSAMOutput' should be cotained in the row names of the input matrix!\n")
    }
    
    module_pc <- apply(hmi[2:nrow(hmi),],1,.prcomp_WGCNA,inputMat,rule)
    module_pc <- t(module_pc)
    rownames(module_pc) <- hmi[2:nrow(hmi),4]
    colnames(module_pc) <- colnames(inputMat)
    
    moduleInfo <-data.frame(moduleName="",featureName="",featureType="",method="",compare="",statistic="",direction=0,pvalue="",fdr="",stringsAsFactors=F)
    mi <- 1
    for(i in c(2:ncol(sampleAnn))){
        
        featureName <- colnames(sampleAnn)[i]
        sampleInfo <- sampleAnn[,i]
        sampleType <- sampleInfo[1]
        sampleInfo <- sampleInfo[2:length(sampleInfo)]
        com <- ""
        sig <- ""
        method <- ""
        if(sampleType=="BIN"){
            re <- .BINTest(module_pc,sampleInfo,BINMethod,fdrmethod,pth)
            method <- BINMethod
            if(!is.null(re)){
                com <- re$dir
                sig <- re$module_sig
            }
        }
        if(sampleType=="CAT"){
            re <- .CATTest(module_pc,sampleInfo,CATMethod,fdrmethod,pth)
            method <- CATMethod
            if(!is.null(re)){
                com <- re$dir
                sig <- re$module_sig
            }

        }
        
        if(sampleType=="CON"){
            sampleInfo <- as.numeric(sampleInfo)
            re <- .CONTest(module_pc,sampleInfo,CONMethod,fdrmethod,pth)
            method <- CONMethod
            if(!is.null(re)){
                com <- ""
                sig <- re
            }
        }
        if(sampleType=="SUR"){
            re <- .SURTest(module_pc,sampleInfo,fdrmethod,pth)
            method <- "Cox model"
            if(!is.null(re)){
                com <- ""
                sig <- re
            }
            
        }
        if(!is.null(nrow(sig))){
            moduleInfo[c(mi:(mi+nrow(sig)-1)),1] <- rownames(sig)
            moduleInfo[c(mi:(mi+nrow(sig)-1)),2] <- rep(featureName,nrow(sig))
            moduleInfo[c(mi:(mi+nrow(sig)-1)),3] <- rep(sampleType,nrow(sig))
            moduleInfo[c(mi:(mi+nrow(sig)-1)),4] <- rep(method,nrow(sig))
            moduleInfo[c(mi:(mi+nrow(sig)-1)),5] <- rep(com,nrow(sig))
            moduleInfo[c(mi:(mi+nrow(sig)-1)),6] <- format(sig[,1],scientific=TRUE, digits=3)
            moduleInfo[c(mi:(mi+nrow(sig)-1)),7] <- sig[,2]
            moduleInfo[c(mi:(mi+nrow(sig)-1)),8] <- format(sig[,3],scientific=TRUE, digits=3)
            moduleInfo[c(mi:(mi+nrow(sig)-1)),9] <- format(sig[,4],scientific=TRUE, digits=3)
            mi <- mi + nrow(sig)
        }
    }
    
    if(moduleInfo[1,1]!=""){
        .createHTML_feature(outputHtmlFile,moduleInfo,pth)
        cat("The assocaitaions between modules and features can be found at ",outputHtmlFile,".html!\n",sep="")
        return (moduleInfo)
    }else{
        cat("There is no assocation between modules and features!\n")
        return(NULL)
    }

}



.prcomp_WGCNA <- function(hmivector,exp,rulef) {
    
    st <- as.numeric(hmivector[5])
    en <- as.numeric(hmivector[6])
    ge <- rulef[st:en,4]
    
    module_exp <- exp[ge,]
    
    col <- rep("black",nrow(module_exp))
    pc <- moduleEigengenes(t(module_exp),col)
    pc <- pc$eigengenes
    
    pc <- as.vector(as.matrix(pc))
    return(pc)
}

.BINTest <- function(module_pc,sampleInfo,method,fdrmethod,pth){
    
    sampleInfo_in <- which(!is.na(sampleInfo))
    sampleInfo <- sampleInfo[sampleInfo_in]
    module_pc_f <- module_pc[,sampleInfo_in]
    
    sampleD <- unique(sampleInfo)
    dir <- paste(sampleD[1],"/",sampleD[2],sep="")
    
    in1 <- which(sampleInfo==sampleD[1])
    in2 <- which(sampleInfo==sampleD[2])
    
    module_sta <- apply(module_pc_f,1,.BINTest_single,in1,in2,method)
    module_sta <- t(module_sta)
    fdrp <- p.adjust(module_sta[,3],method=fdrmethod)
    module_sta <- data.frame(statistic=module_sta[,1],direction=module_sta[,2],pvalue=module_sta[,3],fdr=fdrp,stringsAsFactors=F)
    module_sta_sig <- module_sta[module_sta[,3]<pth,]
    if(nrow(module_sta_sig)>0){
        re <- list(module_sig=module_sta_sig,dir=dir)
        return(re)
    }else{
        return(NULL)
    }
}

.BINTest_single <- function(module_pc_vector,in1,in2,method){
    if(method=="ttest"){
        t <- t.test(module_pc_vector[in1],module_pc_vector[in2])
        statistic <- t$statistic
        direction <- sign(statistic)
        pv <- t$p.value
    }
    if(method=="ranktest"){
        r <- wilcox.test(module_pc_vector[in1],module_pc_vector[in2])
        statistic <- r$statistic
        direction <- ifelse(mean(module_pc_vector[in1])>=mean(module_pc_vector[in2]),1,-1)
        pv <- r$p.value
    }

    re <- c(statistic,direction,pv)
    return(re)

}


.CATTest <- function(module_pc,sampleInfo,method,fdrmethod,pth){
    sampleInfo_in <- which(!is.na(sampleInfo))
    sampleInfo <- sampleInfo[sampleInfo_in]
    module_pc_f <- module_pc[,sampleInfo_in]
    
    sampleD <- unique(sampleInfo)
    if(length(sampleD)<10){
        sampleD <- paste(sampleD,collapse="/")
    }else{
        sampleD <- length(sampleD)
    }
    
    module_sta <- apply(module_pc_f,1,.CATTest_single,sampleInfo,method)
    module_sta <- t(module_sta)
    mediansta <- median(module_sta[,1])
    
    dir <- sign(module_sta[,1]-mediansta)
    fdrp <- p.adjust(module_sta[,2],method=fdrmethod)
    module_sta <- data.frame(statistic=module_sta[,1],direction=dir,pvalue=module_sta[,2],fdr=fdrp,stringsAsFactors=F)
    module_sta_sig <- module_sta[module_sta[,3]<pth,]
    if(nrow(module_sta_sig)>0){
        re <- list(module_sig=module_sta_sig,dir=sampleD)
        return(re)
    }else{
        return(NULL)
    }
}


.CATTest_single <- function(module_pc_vector,group,method){
    
    data <- data.frame(pc=module_pc_vector,group=group)
    
    if(method=="anova"){
        a <- aov(pc ~ group, data=data)
        a <- summary(a)[[1]]
        pv <- a[["Pr(>F)"]][1]
        statistic <- a[["F value"]][1]
    }
    if(method=="kruskal"){
        k <- kruskal.test(pc ~ group,data=data)
        statistic <- k$statistic
        pv <- k$p.value
    }
    
    re <- c(statistic,pv)
    return(re)
    
}

.CONTest <- function(module_pc,sampleInfo,method,fdrmethod,pth){
    sampleInfo_in <- which(!is.na(sampleInfo))
    sampleInfo <- sampleInfo[sampleInfo_in]
    module_pc_f <- module_pc[,sampleInfo_in]
    
    sampleInfo <- as.numeric(sampleInfo)
    
    module_sta <- apply(module_pc_f,1,.CONTest_single,sampleInfo,method)
    module_sta <- t(module_sta)
    fdrp <- p.adjust(module_sta[,3],method=fdrmethod)
    module_sta <- data.frame(statistic=module_sta[,1],direction=module_sta[,2],pvalue=module_sta[,3],fdr=fdrp,stringsAsFactors=F)
    module_sta_sig <- module_sta[module_sta[,3]<pth,]
    if(nrow(module_sta_sig)>0){
        return(module_sta_sig)
    }else{
        return(NULL)
    }
}

.CONTest_single <- function(module_pc_vector,sampleInfo,method){
    c <- cor.test(module_pc_vector,sampleInfo,method=method)
    statistic <- c$estimate
    pv <- c$p.value
    direction <- sign(statistic)
    re <- c(statistic,direction,pv)
    return(re)
}


.SURTest <- function(module_pc,sampleInfo,fdrmethod,pth){
    sampleInfo <- strsplit(sampleInfo,",",fixed=TRUE)
    sampleInfo <- do.call(rbind,sampleInfo)
    time <- sampleInfo[,1]
    event <- sampleInfo[,2]
    time <- suppressWarnings(as.numeric(time))
    event <- suppressWarnings(as.numeric(event))
    time_in <- which(!is.na(time))
    time_f <- time[time_in]
    event_f <- event[time_in]
    module_pc_f <- module_pc[,time_in]
    
    module_sta <- apply(module_pc_f,1,.SURTest_single,time_f,event_f)
    module_sta <- t(module_sta)
    fdrp <- p.adjust(module_sta[,3],method=fdrmethod)
    module_sta <- data.frame(statistic=module_sta[,1],direction=module_sta[,2],pvalue=module_sta[,3],fdr=fdrp,stringsAsFactors=F)
    module_sta_sig <- module_sta[module_sta[,3]<pth,]
    if(nrow(module_sta_sig)>0){
        return(module_sta_sig)
    }else{
        return(NULL)
    }
    
}

.SURTest_single <- function(module_pc_vector,time,event){
    s <- coxph(Surv(time,event) ~ module_pc_vector)
    s <- summary(s)
    statistic <- s$coefficients[1]
    p_ll <- s$logtest[3]
    
    dir <- sign(statistic)
    re <- c(statistic,dir,p_ll)
    return(re)
}


.createHTML_feature <- function(outputHtmlFile,moduleInfo,pth){
    
    outputDir <- dirname(outputHtmlFile)
    htmlFileName <- basename(outputHtmlFile)
    
    target <- HTMLInitFile(outputDir, filename=htmlFileName, BackGroundColor="#BBBBEE")
    HTML(paste("<h2>The associated modules for each feature, p<",pth,"</h2>",sep=""), file=target)
    
    moduleInfo <- moduleInfo[order(moduleInfo[,2],moduleInfo[,8]),]
    
    features <- unique(moduleInfo[,2])
    
    de <- '<table><tr><th width="10%" align="left">Feature</th><th width="10%" align="left">Feature Type</th><th width="10%" align="left">Method</th><th width="20%" align="left">Compare</th><th width="15%" align="left">Module Name</th><th width="10%" align="center">Statistic</th><th width="5%" align="center">Direction</th><th width="10%" align="center">P Value</th><th width="10%" align="center">FDR</th></tr>'
    HTML(de, file=target)
    
    for(i in c(1:length(features))){
        info <- moduleInfo[moduleInfo[,2]==features[i],]
        de <- paste('<tr><td><a name="',info[1,2],'">',info[1,2],"</a></td><td>",info[1,3],"</td><td>",info[1,4],"</td><td>",info[1,5],'</td><td><a href="#',info[1,1],'">',info[1,1],'</a></td><td align="center">',info[1,6],'</td><td align="center">',info[1,7],'</td><td align="center">',info[1,8],'</td><td align="center">',info[1,9],'</td></tr>',sep="")
        HTML(de,file=target)
        if(nrow(info)>1){
            for(j in c(2:nrow(info))){
                de <- paste('<tr><td> </td><td> </td><td> </td><td> </td><td><a href="#',info[j,1],'">',info[j,1],'</a></td><td align="center">',info[j,6],'</td><td align="center">',info[j,7],'</td><td align="center">',info[j,8],'</td><td align="center">',info[j,9],'</td></tr>',sep="")
                HTML(de,file=target)
            }
        }
        
        de <- "<tr><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td></tr>"
        HTML(de,file=target)
        HTML(de,file=target)
    }
    
    HTML("</table><br/><br/>",file=target)
    
    HTML(paste("<h2>The associated features for each module, p<",pth,"</h2>",sep=""), file=target)
    moduleInfo <- moduleInfo[order(moduleInfo[,1],moduleInfo[,8]),]
    modules <- unique(moduleInfo[,1])
    
    de <- '<table><tr><th width="15%" align="left">Module Name</th><th width="10%" align="left">Feature</th><th width="10%" align="left">Feature Type</th><th width="10%" align="left">Method</th><th width="20%" align="left">Compare</th><th width="10%" align="center">Statistic</th><th width="5%" align="center">Direction</th><th width="10%" align="center">P Value</th><th width="10%" align="center">FDR</th></tr>'

    HTML(de, file=target)
    
    for(i in c(1:length(modules))){
        info <- moduleInfo[moduleInfo[,1]==modules[i],]
        de <- paste('<tr><td><a name="',info[1,1],'">',info[1,1],'</a></td><td><a href="#',info[1,2],'">',info[1,2],"</a></td><td>",info[1,3],"</td><td>",info[1,4],"</td><td>",info[1,5],'</td><td align="center">',info[1,6],'</td><td align="center">',info[1,7],'</td><td align="center">',info[1,8],'</td><td align="center">',info[1,9],'</td></tr>',sep="")
        HTML(de,file=target)
        
        if(nrow(info)>1){
            for(j in c(2:nrow(info))){
                de <- paste('<tr><td> </td><td><a href="#',info[j,2],'">',info[j,2],"</a></td><td>",info[j,3],"</td><td>",info[j,4],"</td><td>",info[j,5],'</td><td align="center">',info[j,6],'</td><td align="center">',info[j,7],'</td><td align="center">',info[j,8],'</td><td align="center">',info[j,9],'</td></tr>',sep="")
                HTML(de,file=target)
            }
        }
        
        de <- "<tr><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td></tr>"
        HTML(de,file=target)
        HTML(de,file=target)
        
    }
    HTML("</table><br/><br/>",file=target)
    
    HTMLEndFile()
    
    #browseURL(file.path(outputDir,paste(htmlFileName,".html",sep="")))

}
