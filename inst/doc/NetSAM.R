### R code from vignette source 'NetSAM.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Example
###################################################
library("NetSAM")
inputNetworkDir <- system.file("extdata","exampleNetwork.net",package="NetSAM")
outputFileName <- paste(getwd(),"/NetSAM",sep="")
result <- NetSAM(inputNetwork=inputNetworkDir, outputFileName=outputFileName, outputFormat="nsm", 
edgeType="unweighted", map_to_genesymbol=FALSE, organism="hsapiens", idType="auto", minModule=0.003, 
stepIte=FALSE, maxStep=4, moduleSigMethod="cutoff", modularityThr=0.2, ZRanNum=10, 
PerRanNum=100, ranSig=0.05, edgeThr=(-1), nodeThr=(-1), nThreads=3)



###################################################
### code chunk number 2: Example
###################################################
library("NetSAM")
inputNetworkDir <- system.file("extdata","exampleNetwork.net",package="NetSAM")
outputFileName <- paste(getwd(),"/NetSAM",sep="")
NetAnalyzer(inputNetworkDir,outputFileName,"unweighted")


###################################################
### code chunk number 3: Example
###################################################
library("NetSAM")
inputMatDir <- system.file("extdata","exampleExpressionData_nonsymbol.cct",package="NetSAM")
inputMat <- read.table(inputMatDir,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
mergedData <- mergeDuplicate(id=inputMat[,1],data=inputMat[,2:ncol(inputMat)],collapse_mode="maxSD")


###################################################
### code chunk number 4: Example
###################################################
library("NetSAM")
print("transform ids from a gene list to gene symbols...")
geneListDir <- system.file("extdata","exampleGeneList.txt",package="NetSAM")
geneList <- read.table(geneListDir,header=FALSE,sep="\t",stringsAsFactors=FALSE)
geneList <- as.vector(as.matrix(geneList))
geneList_symbol <- mapToSymbol(inputData=geneList, organism="hsapiens", inputType="genelist",idType="affy_hg_u133_plus_2")
	
print("transform ids in the input network to gene symbols...")
inputNetwork <- system.file("extdata","exampleNetwork_nonsymbol.net",package="NetSAM")
network_symbol <- mapToSymbol(inputData=inputNetwork,organism="hsapiens",inputType="network",idType="entrezgene",edgeType="unweighted")
	
print("transform ids in the input matrix to gene symbols...")
inputMatDir <- system.file("extdata","exampleExpressionData_nonsymbol.cct",package="NetSAM")
matrix_symbol <- mapToSymbol(inputData=inputMatDir,organism="hsapiens",inputType="matrix",idType="affy_hg_u133_plus_2",collapse_mode="maxSD")
	
print("transform ids in the sbt file to gene symbols...")
inputSBTDir <- system.file("extdata","exampleSBT.sbt",package="NetSAM")
sbt_symbol <- mapToSymbol(inputData= inputSBTDir,organism="hsapiens",inputType="sbt",idType="affy_hg_u133_plus_2")
	
print("transform ids in the sct file to gene symbols...")
inputSCTDir <- system.file("extdata","exampleSCT.sct",package="NetSAM")
sct_symbol <- mapToSymbol(inputData= inputSCTDir,organism="hsapiens",inputType="sct",idType="affy_hg_u133_plus_2",collapse_mode="min")


###################################################
### code chunk number 5: Example
###################################################
library("NetSAM")
inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
matNetwork <- MatNet(inputMat=inputMatDir, collapse_mode="maxSD", naPer=0.7, meanPer=0.8, varPer=0.8, 
corrType="spearman", matNetMethod="rank", valueThr=0.6, rankBest=0.003, networkType="signed", 
netFDRMethod="BH", netFDRThr=0.05, idNumThr=(-1), nThreads=3)


###################################################
### code chunk number 6: Example
###################################################
library("NetSAM")
inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
sampleAnnDir <- system.file("extdata","sampleAnnotation.tsi",package="NetSAM")
	
formatedData <- testFileFormat(inputMat=inputMatDir,sampleAnn=sampleAnnDir,collapse_mode="maxSD")


###################################################
### code chunk number 7: Example
###################################################
library("NetSAM")
inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
sampleAnnDir <- system.file("extdata","sampleAnnotation.tsi",package="NetSAM")
data(NetSAMOutput_Example)
outputHtmlFile <- paste(getwd(),"/featureAsso_HTML",sep="")
featureAsso <- featureAssociation(inputMat=inputMatDir, sampleAnn=sampleAnnDir, NetSAMOutput=netsam_output, outputHtmlFile=outputHtmlFile, CONMethod="spearman", CATMethod="kruskal", BINMethod="ranktest", fdrmethod="BH",pth=0.05,collapse_mode="maxSD")


###################################################
### code chunk number 8: Example
###################################################
library("NetSAM")
data(NetSAMOutput_Example)
outputHtmlFile <- paste(getwd(),"/GOAsso_HTML",sep="")
GOAsso <- GOAssociation(NetSAMOutput=netsam_output, outputHtmlFile=outputHtmlFile, organism="hsapiens", fdrmethod="BH", fdrth=0.05, topNum=5)


###################################################
### code chunk number 9: Example
###################################################
library("NetSAM")
inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
sampleAnnDir <- system.file("extdata","sampleAnnotation.tsi",package="NetSAM")
outputFileName <- paste(getwd(),"/MatSAM",sep="")
matModule <- MatSAM(inputMat=inputMatDir, sampleAnn=sampleAnnDir, outputFileName=outputFileName, outputFormat="msm", organism="hsapiens", 
map_to_symbol=FALSE, idType="auto", collapse_mode="maxSD", naPer=0.7, meanPer=0.8, varPer=0.8, corrType="spearman", matNetMethod="rank", 
valueThr=0.6, rankBest=0.003, networkType="signed", netFDRMethod="BH", netFDRThr=0.05, minModule=0.003, stepIte=FALSE, 
maxStep=4, moduleSigMethod="cutoff", modularityThr=0.2, ZRanNum=10, PerRanNum=100, ranSig=0.05, idNumThr=(-1), nThreads=3)


