\name{mapToSymbol}
\alias{mapToSymbol}

\title{
	Map other ids to gene symbols
}

\description{
	 The mapToSymbol function can transform other ids from a gene list, network, matrix, sbt file or sct file to gene symbols.
}

\usage{
mapToSymbol(inputData, organism="hsapiens", inputType="genelist", idType="auto", edgeType="unweighted", collapse_mode="maxSD", is_outputFile=FALSE, outputFileName="", verbose=TRUE)
}

\arguments{
  \item{inputData}{
	mapToSymbol function supports five different types of data:  "genelist", "network", "matrix", "sbt" and "sct".
	For "genelist" type,  \code{inputData} should be a vector containing gene ids. 
	For "network" type, \code{inputData} can be the directory of the input network file including the file name with "net" extension. If \code{edgeType} is "unweighted", each row represents an edge with two node names separated by a tab or space. If \code{edgeType} is "weighted", each row represents an edge with two node names and edge weight separated by a tab or space.  \code{inputNetwork} can also be a data object in R (data object must be igraph, graphNEL, matrix or data.frame class). 
	For "matrix" type, \code{inputData} should be a directory containing a file name with extension "cct" or "cbt"  or a matrix or data.frame object in R. The first column and first row of the "cct" or "cbt" file should be the row and column names, respectively and other parts are the numeric values. The detail information of "cct" or "cbt" format can be found in the manual of NetGestalt (www.netgestalt.org). A matrix or data.frame object should have row and column names and only contain numeric or integer values.
	For "sbt" type, \code{inputData} should be a directory containing a file name with extension "sbt" . The first column of a "sbt" file is the track names, the second one is the descriptions and others are the ids contained in the track. A "sbt" file can contain multiple tracks. The detail information of "sbt" format can be found in the manual of NetGestalt (www.netgestalt.org). 
	For "sct" type, \code{inputData} should be a directory containing a file name with extension "sct" . The first column of a "sct" file is id names, the first row is the column names and others are the numeric or integer values. The detail information of "sct" format can be found in the manual of NetGestalt (www.netgestalt.org). 
  }
  \item{organism}{
	The organism of the input data. Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana. The default is "hsapiens". 
  }
  \item{inputType}{
  	The type of the input data. see detail information in \code{inputData}.
  }
   \item{idType}{
  	The id type of the ids in the input data. MatSAM will use BiomaRt package to transform the input ids to gene symbols based on \code{idType}. The users can also set \code{idType} as "auto" that means MatSAM will automatically search the id type of the input data. However, this may take 10 minutes based on the users' internet speed. The default is "auto".
  }
  \item{edgeType}{
	The type of the input network: "weighted" or "unweighted". 
  }
 \item{collapse_mode}{
	The method to collapse duplicate ids. "mean", "median", "maxSD", "maxIQR", "max" and "min" represent the mean, median, max standard deviation, max interquartile range, maximum and minimum of values for ids in each sample. The default is "maxSD". For SCT file, we suggest to use "max" or "min" to collapse duplicate ids in the statistic data.
  }
  \item{is_outputFile}{
  	If \code{is_outputFile} is TRUE, the function will output the transformed data to a file. The default is FALSE.
  }
    \item{outputFileName}{
  	The output file name.
  }
  \item{verbose}{
	Report the extra information on progress. The default is TRUE.
  }
}

\author{
	Jing Wang
}

\examples{
	
	###transform ids from a gene list to gene symbols###
	geneListDir <- system.file("extdata","exampleGeneList.txt",package="NetSAM")
	geneList <- read.table(geneListDir,header=FALSE,sep="\t",stringsAsFactors=FALSE)
	geneList <- as.vector(as.matrix(geneList))
	geneList_symbol <- mapToSymbol(inputData=geneList, organism="hsapiens", inputType="genelist",idType="affy_hg_u133_plus_2")
	
	###transform ids in the input network to gene symbols###
	inputNetwork <- system.file("extdata","exampleNetwork_nonsymbol.net",package="NetSAM")
	network_symbol <- mapToSymbol(inputData=inputNetwork,organism="hsapiens",inputType="network",idType="entrezgene",edgeType="unweighted")
	
	###transform ids in the input matrix to gene symbols###
	inputMatDir <- system.file("extdata","exampleExpressionData_nonsymbol.cct",package="NetSAM")
	matrix_symbol <- mapToSymbol(inputData=inputMatDir,organism="hsapiens",inputType="matrix",idType="affy_hg_u133_plus_2",collapse_mode="maxSD")
	
	###transform ids in the sbt file to gene symbols###
	inputSBTDir <- system.file("extdata","exampleSBT.sbt",package="NetSAM")
	sbt_symbol <- mapToSymbol(inputData= inputSBTDir,organism="hsapiens",inputType="sbt",idType="affy_hg_u133_plus_2")
	
	###transform ids in the sct file to gene symbols###
	inputSCTDir <- system.file("extdata","exampleSCT.sct",package="NetSAM")
	sct_symbol <- mapToSymbol(inputData= inputSCTDir,organism="hsapiens",inputType="sct",idType="affy_hg_u133_plus_2",collapse_mode="min")

}

\keyword{ methods }
