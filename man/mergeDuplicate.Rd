\name{mergeDuplicate}
\alias{mergeDuplicate}

\title{
	Merge the duplicate Ids in the matrix data
}

\description{
	 The mergeDuplicate function will merge the duplicate Ids in the matrix data and return the matrix with unique Ids. This function can also used to merge the duplicate mapped Ids when transforming the Ids of data matrix to other Ids.

}

\usage{
	mergeDuplicate(id, data, collapse_mode="maxSD")
}

\arguments{
  \item{id}{
	Duplicate Ids that should be a vector object in R.
  }
  \item{data}{
  	the corresponding data matrix that has the same number of rows with \code{id} and should be a matrix or data.frame object in R.
  }
 \item{collapse_mode}{
	The method to collapse duplicate ids. "mean", "median", "maxSD", "maxIQR", "max" and "min" represent the mean, median, max standard deviation, max interquartile range, maximum and minimum of values for ids in each sample. The default is "maxSD".
  }

}

\value{
	The function will return the data matrix with unique Ids.
}

\author{
	Jing Wang
}

\examples{
	inputMatDir <- system.file("extdata","exampleExpressionData_nonsymbol.cct",package="NetSAM")
	inputMat <- read.table(inputMatDir,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
	mergedData <- mergeDuplicate(id=inputMat[,1],data=inputMat[,2:ncol(inputMat)],collapse_mode="maxSD")
}

\keyword{ methods }
