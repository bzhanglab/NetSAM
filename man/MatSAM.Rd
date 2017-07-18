\name{MatSAM}
\alias{MatSAM}

\title{
	Correlation network construction, seriation and modularization from a matrix
}

\description{
	 The MatSAM function first uses MatNet function to identify the correlation network and then uses NetSAM function to identify the module and optimize the one-dimensional ordering of the nodes in each module. 
}

\usage{
MatSAM(inputMat, sampleAnn=NULL, outputFileName, outputFormat="msm", organism="hsapiens", map_to_symbol=FALSE, idType="auto", collapse_mode="maxSD", naPer=0.7, meanPer=0.8, varPer=0.8, corrType="spearman", matNetMethod="rank", valueThr=0.5, rankBest=0.003, networkType="signed", netFDRMethod="BH", netFDRThr=0.05, minModule=0.003, stepIte=FALSE, maxStep=4, moduleSigMethod="cutoff", modularityThr=0.2, ZRanNum=10, PerRanNum=100, ranSig=0.05,  idNumThr=(-1), nThreads=3)
}

\arguments{
  \item{inputMat}{
	\code{inputMat} should contain a file name with extension "cct" or "cbt" or a matrix or data.frame object in R. The first column and first row of the "cct" or "cbt" file should be the row and column names, respectively and other parts are the numeric values. The detail information of "cct" or "cbt" format can be found in the manual of NetGestalt (www.netgestalt.org). A matrix or data.frame object should have row and column names and only contain numeric or integer values.
  }
  \item{sampleAnn}{
	\code{sampleAnn} should contain a file name with  "tsi" extension (the detail information of "tsi" format can be found in the manual of NetGestalt (www.netgestalt.org)) or a data.frame object in R. If the data does not have sample annotation, this argument can be ignored. The first row of the data is the name  of sample features. The second row is the type of each feature. The third row is the category of each feature. If there is no category information for the features, the sample information will start from the third row . The first column is the sample name. 
  }
  \item{outputFileName}{
	Output file name. The file name extension is "msm" which can be uploaded to the NetGestalt directly.
  }
  \item{outputFormat}{
	The format of the output file. "msm" format can be used as an input in NetGestalt; "gmt" format can be used to do other network analysis (e.g. as an input in GSEA (Gene Set Enrichment Analysis) to do module enrichment analysis); "multiple" represents the MatSAM function will output five files: ruler file containing gene order information, hmi file containing module information, net file containing correlation network information, cct file containing the filtered data matrix, and tsi file containing the sample annotation with standardized format; and "none" represents the function will not output any file.
  }
  \item{organism}{
	The organism of the input data. Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana. The default is "hsapiens". 
  }
  \item{map_to_symbol}{
	If \code{map_to_symbol} is TRUE, the function will first change the input ids to gene symbols and collapse multiple ids with the same gene symbol based on the \code{collapse_mode} method before identifying correlation network. The default is FALSE.
  }
    \item{idType}{
  	The id type of the ids in the input matrix. MatSAM will use BiomaRt package to transform the input ids to gene symbols based on \code{idType}. The users can also set \code{idType} as "auto" that means MatSAM will automatically search the id type of the input data. However, this may take 10 minutes based on the users' internet speed. The default is "auto".
  }
 \item{collapse_mode}{
	The method to collapse duplicate ids. "mean", "median", "maxSD", "maxIQR", "max" and "min" represent the mean, median, max standard deviation, max interquartile range, maximum and minimum of values for ids in each sample. The default is "maxSD".
  }
  \item{naPer}{
	To remove ids with missing values in most of samples, the function calculates the percentage of missing values in all samples for each id and removes ids with over \code{naPer} missing values in all samples. The default \code{naPer} is 0.7.
 }
  \item{meanPer}{
	To remove ids with low values, the function calculates the mean of values for a id in all samples and remains top \code{meanPer} ids based on the mean. The default \code{meanPer} is 0.8.
  }
  \item{varPer}{
	Based on the remained ids filtered by \code{meanPer}, the function can also remove less variable ids by calculating the standard deviation of values for a id in all samples and remaining top \code{varPer} ids based on the standard deviation. The default \code{varPer} is 0.8.
  }
  \item{corrType}{
	A character string indicating which correlation coefficient is to be computed for each pair of ids. The function supports "spearman" (default) or "pearson" method.
  }
  \item{matNetMethod}{
	MatNet function supports three methods to construct correlation network: "value", "rank" and "directed". 1. "value" method: the correlation network only remains id pairs with correlations over cutoff threshold \code{valueThr}; 2. "rank" method: for each id A, the function first selects ids that significantly correlate with id A and then extracts a set of ids (the number of ids is calculated based on \code{rankBest}) that are most similar to id A from the significant set. Then, for each id B in the set, the function also extracts the same number of ids that are significant correlated and most similar to id B. If id A is in the set of id B, the edge between id A and id B will be remained. Combining all remained edges can construct a correlation network; 3. "directed" method: the function will only remain the best significant id for each id as the edge.Combining all edges can construct a directed correlation network.
  }
  \item{valueThr}{
	Correlation cutoff threshold for "value" method. The default is 0.5.
  }
  \item{rankBest}{
	The percentage of ids that are most similar to one id for "rank" method. The default is 0.003 which means the "rank" method will select top 30 most similar ids for each id if the number of ids in the matrix is 10,000.
  }
  \item{networkType}{
	If \code{networkType} is "unsigned", the correlation of all pairs of ids will be changed to absolute values. The default is "signed".
  }
  \item{netFDRMethod}{
	p value adjustment methods for "rank" and "directed" methods. The default is "BH".
  }
  \item{netFDRThr}{
	fdr threshold for identifying significant pairs for "rank" and "directed" methods. The default is 0.05
  }
\item{minModule}{
	The minimum percentage of nodes in a module. The minimum size of a module is calculated by multiplying \code{minModule} by the number of nodes in the whole network. If the size of a module identified by the function is less than the minimum size, the module will not be further partitioned into sub-modules. The default is 0.003 which means the minimum module size is 30 if there are 10,000 nodes in the whole network. If the minimum module size is less than 5, the minimum module size will be set as 5. The \code{minModule} should be less than 0.2.
  }
  \item{stepIte}{
	Because NetSAM uses random walk distance-based hierarchical clustering to reveal the hierarchical organization of an input network, it requires a specified length of the random walks. If \code{stepIte} is TRUE, the function will test a range of lengths ranging from 2 to \code{maxStep} to get the optimal length. Otherwise, the function will directly use \code{maxStep} as the length of the random walks. The default \code{maxStep} is 4. Because optimizing the length of the random walks will take a long time, if the network is too big (e.g. the number of edges is over 200,000), we suggest to set \code{stepIte} as FALSE.
  }
  \item{maxStep}{
	The length or max length of the random walks.
  }
  \item{moduleSigMethod}{
	To test whether a network under consideration has a non-random internal modular organization, the function provides three options: "cutoff", "zscore" and "permutation". "cutoff" means if the modularity score of the network is above a specified cutoff value, the network will be considered to have internal organization and will be further partitioned. For "zscore" and "permutation", the function will first generate a set of random modularity scores. Based on a unweighted network, the function uses the edge switching method to generate a given number of random networks with the same number of nodes and an identical degree sequence and calculates the modularity scores for these random networks. Based on a weighted network, the function shuffles the weights of all edges and calculate the modularity scores for network with random weights. Then, "zscore" method will transform the real modularity score to a z score based on the random modularity scores and then transform the z score to a p value assuming a standard normal distribution. The "permutation" method will compare the real modularity score with the random ones to calculate a p value. Finally, under a specified significance level, the function determines whether the network can be further partitioned. The default is "cutoff".
  }
  \item{modularityThr}{
	Threshold of modularity score for the  "cutoff" method. The default is 0.2
  }
  \item{ZRanNum}{
	The number of random networks that will be generated for the "zscore" calculation. The default is 10.
  }
  \item{PerRanNum}{
	The number of random networks that will be generated for the "permutation" p value calculation. The default is 100.
  }
  \item{ranSig}{
	The significance level for determining whether a network has non-random internal modular organization for the "zscore" or "permutation" methods. The default is 0.05.
  }
  \item{idNumThr}{
  	If the matrix contains too many ids, it will take a long time and use a lot of memory to identify the modules. Thus, the function provides the option to set the threshold of number of ids for further analysis. After filtering by meanPer and varPer, if the number of ids is still larger than \code{idNumThr}, the function will select top \code{idNumThr} ids with the largest variance. The default is -1, which means there is no limitation for the matrix.
  }
  \item{nThreads}{
	MatSAM function supports parallel computing based on multiple cores. The default is 3.
  }
}

\value{
	Including a "msm" file, the function will output a list object containing module information, gene order information, correlation network and filtered matrix based on the ids in the network. The function will also output two HTML files that contain the significant associations between sample features and modules and associated GO terms for the modules.
}

\note{
After identifying the modules, the MatSAM function will identify the associations between sample features and modules using the featureAssociation function or the associated GO terms for the modules using the GOAssociation function. For the featureAssociation function, MatSAM only uses the default parameters. For the GOAssociation function, MatSAM sets "outputType" as "top" and "topNum" as 1. The users can use the list object returned by MatSAM as the input of the function featureAssociation and GOAssociation to perform some further analysis based on the different parameters.
}


\author{
	Jing Wang
}

\examples{
	inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
	cat(inputMatDir)
	sampleAnnDir <- system.file("extdata","sampleAnnotation.tsi",package="NetSAM")
	cat(sampleAnnDir)
	outputFileName <- paste(getwd(),"/MatSAM",sep="")
	matModule <- MatSAM(inputMat=inputMatDir, sampleAnn=sampleAnnDir, outputFileName=outputFileName, outputFormat="msm", organism="hsapiens", map_to_symbol=FALSE, idType="auto", collapse_mode="maxSD", naPer=0.7, meanPer=0.8, varPer=0.8, corrType="spearman", matNetMethod="rank", valueThr=0.6, rankBest=0.003, networkType="signed", netFDRMethod="BH", netFDRThr=0.05, minModule=0.003, stepIte=FALSE, maxStep=4, moduleSigMethod="cutoff", modularityThr=0.2, ZRanNum=10, PerRanNum=100, ranSig=0.05, idNumThr=(-1), nThreads=3)
}

\seealso{
	\code{\link{MatNet}}
	\code{\link{NetSAM}}
}

\keyword{ methods }
