\name{MatNet}
\alias{MatNet}

\title{
	Construction of correlation network from a matrix
}

\description{
	 The MatNet function can use one of three different methods to construct correlation network based on the input data matrix. The output correlation network can be used as an input of NetSAM function to identify hierarchical modules.
}

\usage{
MatNet(inputMat, collapse_mode="maxSD", naPer=0.7, meanPer=0.8, varPer=0.8, corrType="spearman", matNetMethod="rank", valueThr=0.5, rankBest=0.003, networkType="signed", netFDRMethod="BH", netFDRThr=0.05, idNumThr=(-1),nThreads=3)
}

\arguments{
 \item{inputMat}{
	\code{inputMat} should contain a file name with extension "cct" or "cbt" or a matrix or data.frame object in R. The first column and first row of the "cct" or "cbt" file should be the row and column names, respectively and other parts are the numeric values. The detail information of "cct" or "cbt" format can be found in the manual of NetGestalt (www.netgestalt.org). A matrix or data.frame object should have row and column names and only contain numeric or integer values.
  }
  \item{collapse_mode}{
	If the input matrix data contains the duplicate ids, the function will collapse duplicate ids based on the \code{collapse_mode}. "mean", "median", "maxSD" and "maxIQR" represent the mean, median, max standard deviation or max interquartile range of id values in each sample. The default is "maxSD".
  }
   \item{naPer}{
	To remove ids with missing values in most of samples, the function calculates the percentage of missing values in all samples for each id and removes ids with over \code{naPer} missing values in all samples. The default \code{naPer} is 0.7.
 }
 \item{meanPer}{
	To remove ids with low values, the function calculates the mean of values for each id in all samples and remains top \code{meanPer} ids based on the mean. The default \code{meanPer} is 0.8.
  }
  \item{varPer}{
	Based on the remained ids filtered by \code{meanPer}, the function can also remove less variable ids by calculating the standard deviation of values for each id in all samples and remaining top \code{varPer} ids based on the standard deviation. The default \code{varPer} is 0.8.
  }
  \item{corrType}{
	The method to calculate correlation coefficient for each pair of ids. The function supports "spearman" (default) or "pearson" method.
  }
  \item{matNetMethod}{
	MatNet function supports three methods to construct correlation network: "value", "rank" and "directed". 1. "value" method: the correlation network only remains id pairs with correlations over cutoff threshold \code{valueThr}; 2. "rank" method: for each id A, the function first selects ids that significantly correlate with id A and then extracts a set of candidate neighbors (the number of ids is calculated based on \code{rankBest}) from the significant set that are most similar to id A. Then, for each id B in the candidate neighbors of id A, the function also extracts the same number of ids that are significant correlated and most similar to id B. If id A is also the candidate neighbors of id B, there will be an edge between id A and id B. Combining all edges can construct a correlation network; 3. "directed" method: the function will only remain the best significant id for each id as the edge.Combining all edges can construct a directed correlation network.
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
  \item{idNumThr}{
  	If the matrix contains too many ids, it will take a long time and use a lot of memory to identify the modules. Thus, the function provides the option to set the threshold of number of ids for further analysis. After filtering by meanPer and varPer, if the number of ids is still larger than \code{idNumThr}, the function will select top \code{idNumThr} ids with the largest variance. The default is -1, which means there is no limitation for the matrix.
  }
  \item{nThreads}{
	MatNet function supports parallel computing based on multiple cores. The default is 3.
  }
}

\note{
For data with missing values, the function will take longer time to calculate correlation between each pair of ids than data without missing value.
}

\author{
	Jing Wang
}

\examples{
	inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
	matNetwork <- MatNet(inputMat=inputMatDir, collapse_mode="maxSD", naPer=0.7, meanPer=0.8, varPer=0.8, corrType="spearman", matNetMethod="rank", valueThr=0.6, rankBest=0.003, networkType="signed", netFDRMethod="BH", netFDRThr=0.05, idNumThr=(-1),nThreads=3)
}

\keyword{ methods }
