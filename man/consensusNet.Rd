\name{consensusNet}
\alias{consensusNet}

\title{
	Construction of a consensus coexpression network
}

\description{
  To increase robustness against errors in data, the consensusNet function uses a bootstrapping procedure to construct a coexpression network.
}

\usage{
  consensusNet(data, organism="hsapiens",bootstrapNum=100, naPer=0.5, meanPer=0.8,varPer=0.8,method="rank_unsig",value=3/1000,pth=1e-6, nMatNet=2, nThreads=4)
}

\arguments{
 \item{data}{
	\code{data} should contain a file name with extension "cct" or "cbt" or a matrix or data.frame object in R. The first column and first row of the "cct" or "cbt" file should be the row and column names, respectively and other parts are the numeric values. The detail information of "cct" or "cbt" format can be found in the manual of NetGestalt (www.netgestalt.org). A matrix or data.frame object should have row and column names and only contain numeric or integer values.
  }
 \item{organism}{
   	The organism of the input data. Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, 
celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana. The default is "hsapiens".
 }
 \item{bootstrapNum}{
   Number of bootstrap data sets generated. Default is 100.
 }
   \item{naPer}{
	To remove ids with missing values in most of samples, the function calculates the percentage of missing values in all samples for each id and removes ids with over \code{naPer} missing values in all samples. The default \code{naPer} is 0.5.
 }
 \item{meanPer}{
	To remove ids with low values, the function calculates the mean of values for each id in all samples and remains top \code{meanPer} ids based on the mean. 
The default \code{meanPer} is 0.8.
  }
  \item{varPer}{
	Based on the remaining ids filtered by \code{meanPer}, the function can also remove less variable ids by calculating the standard deviation of values for each id in all
 samples and remaining top \code{varPer} ids based on the standard deviation. The default \code{varPer} is 0.8.
  }
  \item{method}{
  Method used for constructing correlation network with \code{MatNet}. Currently supports "rank", "value" and "rank_unsig". Default is "rank_unsig".
  }
  \item{value}{
  The corresponding value set for \code{method}. Default is 0.003.
  }
  \item{pth}{
  p-value threshold for including an edge. Default is 1.0e-6.
  }
  \item{nMatNet}{
   The number of concurrent running MatNet processes, default is 2.
  }
  \item{nThreads}{
	 consensusNet function supports parallel computing based on multiple cores. The default is 4.
  }
}

\author{
	Jing Wang
}

\examples{
  inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
  data <- read.table(inputMatDir, header=TRUE, row.names=1, stringsAsFactors=FALSE)
  net <- consensusNet(data, organism="hsapiens",bootstrapNum=10, naPer=0.5, meanPer=0.8,varPer=0.8,method="rank_unsig",value=3/1000,pth=1e-6, nMatNet=2, nThreads=4)
}

\keyword{ methods }
