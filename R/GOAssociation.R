GOAssociation <-
function(NetSAMOutput,outputHtmlFile,organism="hsapiens",outputType="significant",fdrmethod="BH",fdrth=0.05,topNum=5){
    
    #require("GO.db") || stop("package GO.db is required")
    
    if(is.null(NetSAMOutput$rulfile) || is.null(NetSAMOutput$hmifile)){
        stop("NetSAMOutput should be the output from NetSAM function, which contain rulfile, hmifile and network!\n")
    }
    
    if(length(which(outputType %in% c("significant","top")))==0){
        stop("The input 'outputType' is invalid! Please select a method from 'significant' and 'top'!\n")
    }
    
    organisms <- c("hsapiens","mmusculus","rnorvegicus","drerio","celegans","scerevisiae","cfamiliaris","dmelanogaster","athaliana")
    names(organisms) <- c("Hs","Mm","Rn","Dr","Ce","Sc","Cf","Dm","At")
    
    if(length(which(organisms==organism))==0){
        stop("Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana!")
    }

    or <- names(organisms)[which(organisms==organism)]
    
    #if(organism=="scerevisiae"){
    #    require("org.Sc.sgd.db") || stop("package org.Sc.sgd.db is required!!\n")
    #}else{
    #    if(organism=="athaliana"){
    #        require("org.At.tair.db") || stop("package org.At.tair.db is required!!\n")
    #    }else{
    #        require(paste("org.",or,".eg.db",sep=""),character.only = TRUE) || stop(paste("package org.", or, ".eg.db is required!!\n", sep=""))
    #    }

    #}
    
    if(length(which(fdrmethod %in% c("holm","hochberg","hommel","bonferroni","BH","BY","none")))==0){
        stop("The input 'fdrmethod' is invalid! Please select an option from 'holm','hochberg', 'hommel', 'bonferroni', 'BH', 'BY' and 'none'!\n")
    }
    
    if(organism=="scerevisiae"){
        .sql_bp <-  "select distinct t1.gene_name,t2.go_id from sgd as t1 inner join go_bp_all as t2 on t1._id=t2._id"
        .sql_cc <-  "select distinct t1.gene_name,t2.go_id from sgd as t1 inner join go_cc_all as t2 on t1._id=t2._id"
        .sql_mf <-  "select distinct t1.gene_name,t2.go_id from sgd as t1 inner join go_mf_all as t2 on t1._id=t2._id"
    }else{
        .sql_bp <-  "select distinct t1.symbol,t2.go_id from gene_info as t1 inner join go_bp_all as t2 on t1._id=t2._id"
        .sql_cc <-  "select distinct t1.symbol,t2.go_id from gene_info as t1 inner join go_cc_all as t2 on t1._id=t2._id"
        .sql_mf <-  "select distinct t1.symbol,t2.go_id from gene_info as t1 inner join go_mf_all as t2 on t1._id=t2._id"
    }
    if(organism=="scerevisiae"){
        conn <- org.Sc.sgd.db::org.Sc.sgd_dbconn()
    }else{
        if(organism=="athaliana"){
            conn <- org.At.tair.db::org.At.tair_dbconn()
        }else{
            conn <- eval(parse(text=paste("org.",or,".eg.db::","org.",or,".eg_dbconn()", sep = "")))
        }
    }
    generalAnn_BP <- dbGetQuery(conn, .sql_bp)
    generalAnn_BP <- generalAnn_BP[!is.na(generalAnn_BP[,1]),]
    bp_num <- tapply(generalAnn_BP[,1],generalAnn_BP[,2],length)
    bp_num <- bp_num[bp_num>=10 & bp_num<=2000]
    generalAnn_BP <- generalAnn_BP[generalAnn_BP[,2] %in% names(bp_num),]
    
    generalAnn_CC <- dbGetQuery(conn, .sql_cc)
    generalAnn_CC <- generalAnn_CC[!is.na(generalAnn_CC[,1]),]
    cc_num <- tapply(generalAnn_CC[,1],generalAnn_CC[,2],length)
    cc_num <- cc_num[cc_num>=10 & cc_num<=2000]
    generalAnn_CC <- generalAnn_CC[generalAnn_CC[,2] %in% names(cc_num),]
    
    generalAnn_MF <- dbGetQuery(conn, .sql_mf)
    generalAnn_MF <- generalAnn_MF[!is.na(generalAnn_MF[,1]),]
    mf_num <- tapply(generalAnn_MF[,1],generalAnn_MF[,2],length)
    mf_num <- mf_num[mf_num>=10 & mf_num<=2000]
    generalAnn_MF <- generalAnn_MF[generalAnn_MF[,2] %in% names(mf_num),]
    
    
    .sql <- "select distinct go_id goid,term name from go_term";
    conn <- get("GO_dbconn")()
    allTermName <- dbGetQuery(conn,.sql)



    rul <- NetSAMOutput$rulfile
    hmi <- NetSAMOutput$hmifile
    
    refGenes <- rul[,4]
    
    annRef_BP <- generalAnn_BP[generalAnn_BP[,1] %in% refGenes,]
    annRef_CC<- generalAnn_CC[generalAnn_CC[,1] %in% refGenes,]
    annRef_MF <- generalAnn_MF[generalAnn_MF[,1] %in% refGenes,]
    
    moduleEnrich <- data.frame(moduleName="",Ontology="",GOID="",GOName="",pvalue=0,fdr=0,stringsAsFactors=F)
    mi <- 1
    for(i in c(2:nrow(hmi))){
        mN <- hmi[i,4]
        s <- hmi[i,5]
        e <- hmi[i,6]
        
        g <- rul[s:e,4]
        
        annInterest_BP <- generalAnn_BP[generalAnn_BP[,1] %in% g,]
        annInterest_CC <- generalAnn_CC[generalAnn_CC[,1] %in% g,]
        annInterest_MF <- generalAnn_MF[generalAnn_MF[,1] %in% g,]
        
        t <- c()
        if(nrow(annInterest_BP)>0){
            termInfo <- .enrichmentFunction(annRef_BP, annInterest_BP,allTermName,fdrmethod)
            if(outputType=="significant"){
                termInfo <- termInfo[termInfo[,6]<fdrth,]
                if(nrow(termInfo)>0){
                    termInfo$ontology <- "Biological Process"
                    t <- rbind(t,termInfo)
                }
            }else{
                termInfo <- termInfo[order(termInfo[,5]),]
                termInfo <- termInfo[1:topNum,]
                termInfo$ontology <- "Biological Process"
                
                t <- rbind(t,termInfo)
            }
        }
        if(nrow(annInterest_CC)>0){
            termInfo <- .enrichmentFunction(annRef_CC, annInterest_CC,allTermName,fdrmethod)
            if(outputType=="significant"){
                termInfo <- termInfo[termInfo[,6]<fdrth,]
                if(nrow(termInfo)>0){
                    termInfo$ontology <- "Cellular Component"
                    t <- rbind(t,termInfo)
                }
            }else{
                termInfo <- termInfo[order(termInfo[,5]),]
                termInfo <- termInfo[1:topNum,]
                termInfo$ontology <- "Cellular Component"
                t <- rbind(t,termInfo)
            }
        }
        if(nrow(annInterest_MF)>0){
            termInfo <- .enrichmentFunction(annRef_MF, annInterest_MF,allTermName,fdrmethod)
            if(outputType=="significant"){
                termInfo <- termInfo[termInfo[,6]<fdrth,]
                if(nrow(termInfo)>0){
                    termInfo$ontology <- "Molecular Function"
                    t <- rbind(t,termInfo)
                }
            }else{
                termInfo <- termInfo[order(termInfo[,5]),]
                termInfo <- termInfo[1:topNum,]
                termInfo$ontology <- "Molecular Function"
                t <- rbind(t,termInfo)
            }
        }
        
        if(length(t)>0){
            moduleEnrich[mi:(mi+nrow(t)-1),1] <- mN
            moduleEnrich[mi:(mi+nrow(t)-1),2] <- t[,7]
            moduleEnrich[mi:(mi+nrow(t)-1),3] <- t[,1]
            moduleEnrich[mi:(mi+nrow(t)-1),4] <- t[,2]
            moduleEnrich[mi:(mi+nrow(t)-1),5] <- format(t[,5],scientific=TRUE,digits=3)
            moduleEnrich[mi:(mi+nrow(t)-1),6] <- format(t[,6],scientific=TRUE,digits=3)
            mi <- mi + nrow(t)
        }

    }

    if(moduleEnrich[1,1]!=""){
        .createHTML_GO(outputHtmlFile,moduleEnrich,outputType,fdrmethod,fdrth,topNum)
        cat("The associated GO terms for each module can be found at ",outputHtmlFile,".html!\n",sep="")
        return(moduleEnrich)
    }else{
        cat("There is no associated GO term for each module based on the input parameters!\n")
        return(NULL)
    }

}


.enrichmentFunction <- function(annRef, annInterest, allTermName, fdrmethod)
{
    allRefnum <- length(unique(annRef[,1]))
    allInterestnum <- length(unique(annInterest[,1]))
    
    allAnnterm <- unique(annRef[,2])
    allAnntermL <- length(allAnnterm)
    
    
    refTermCount <- tapply(annRef[,1],annRef[,2],length)
    
    refTermName <- allTermName[allTermName[,1] %in% unique(annRef[,2]),]
    rownames(refTermName) <- refTermName[,1]
    refTermName <- refTermName[names(refTermCount),]
    
    refTermCount <- data.frame(goid=names(refTermCount),name=refTermName[,2],refnum=refTermCount,stringsAsFactors=F)
    refTermCount <- refTermCount[order(refTermCount[,1]),]
    interestTermCount <- tapply(annInterest[,1],annInterest[,2],length)
    interestTermCount <- data.frame(goid=names(interestTermCount),interestnum=interestTermCount,stringsAsFactors=F)
    interestTermCount <- interestTermCount[order(interestTermCount[,1]),]
    
    ref_interest_TermCount <- refTermCount
    
    ref_interest_TermCount$interestnum = array(0,dim=c(length(ref_interest_TermCount$goid),1))
    ref_interest_TermCount[ref_interest_TermCount$goid %in% interestTermCount[,1],4]=interestTermCount$interestnum
    
    
    n <- nrow(ref_interest_TermCount)
    pv <- array(0,dim=c(n,1))
    for (i in c(1:n)){
        p <- 1-phyper(ref_interest_TermCount[i,4]-1,allInterestnum,allRefnum-allInterestnum,ref_interest_TermCount[i,3],lower.tail = TRUE,log.p= FALSE)
        pv[i,1] <- p
    }
    ref_interest_TermCount$pvalue <- pv
    adp <- p.adjust(pv,method=fdrmethod)
    ref_interest_TermCount$FDR <- adp
    return(ref_interest_TermCount)
}





.createHTML_GO <- function(outputHtmlFile,moduleEnrich,outputType,fdrmethod,fdrth,topNum){
    
    outputDir <- dirname(outputHtmlFile)
    htmlFileName <- basename(outputHtmlFile)
    
    target <- HTMLInitFile(outputDir, filename=htmlFileName, BackGroundColor="#BBBBEE")
    HTML("<h2>The associated GO terms for each module", file=target)
    if(outputType=="significant"){
        HTML(paste("<h2>The GO terms for each Ontology were identified under FDR ",fdrth," based on ",fdrmethod,"</h2>",sep=""), file=target)
    }else{
        HTML(paste("<h2>Top ",topNum," GO terms for each Ontology were selected as the related GO terms</h2>",sep=""), file=target)
    }
    
    moduleEnrich <- moduleEnrich[order(moduleEnrich[,1],moduleEnrich[,2],moduleEnrich[,5]),]
    
    modules <- unique(moduleEnrich[,1])
    
    de <- '<table><tr><th width="15%" align="left">Module Name</th><th width="15%" align="left">Ontology</th><th width="15%" align="left">GO ID</th><th width="35%" align="left">GO Name</th><th width="10%" align="center">P Value</th><th width="10%" align="center">FDR</th></tr>'

    HTML(de, file=target)
    
    for(i in c(1:length(modules))){
        info <- moduleEnrich[moduleEnrich[,1]==modules[i],]
        gohref <- paste("http://amigo.geneontology.org/amigo/term/",info[i,3],sep="")
        de <- paste('<tr><td>',info[1,1],'</td><td>',info[1,2],'</td><td><a href="',gohref,'" target="_blank">',info[1,3],"</a></td><td>",info[1,4],'</td><td align="center">',info[1,5],'</td><td align="center">',info[1,6],"</td></tr>",sep="")
        HTML(de,file=target)
        
        if(nrow(info)>1){
            for(j in c(2:nrow(info))){
                 de <- paste('<tr><td> </td><td>',info[j,2],'</td><td><a href="',gohref,'" target="_blank">',info[j,3],"</a></td><td>",info[j,4],'</td><td align="center">',info[j,5],'</td><td align="center">',info[j,6],"</td></tr>",sep="")
                 HTML(de,file=target)
            }
        }
        
        de <- "<tr><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td><td> </td></tr>"
        HTML(de,file=target)
        HTML(de,file=target)
        
    }
    HTML("</table><br/><br/>",file=target)
    
    HTMLEndFile()
    
    #browseURL(file.path(outputDir,paste(htmlFileName,".html",sep="")))

}
