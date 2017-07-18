\name{testFileFormat}
\alias{testFileFormat}

\title{
	Test whether the data matrix and the annotation have a correct format
}

\description{
	 The testFileFormat function will test the format of the input data matrix and annotation data and return the standardized data matrix and sample annotation data.
}

\usage{
	testFileFormat(inputMat=NULL,sampleAnn=NULL,collapse_mode="maxSD")
}

\arguments{

   \item{inputMat}{
	\code{inputMat} should contain a file name with extension "cct" or "cbt" or a matrix or data.frame object in R. The first column and first row of the "cct" or "cbt" file should be the row and column names, respectively and other parts are the numeric values. The detail information of "cct" or "cbt" format can be found in the manual of NetGestalt (www.netgestalt.org). A matrix or data.frame object should have row and column names and only contain numeric or integer values.
  }
  \item{sampleAnn}{
	\code{sampleAnn} is a file name or a data.frame object in R.
  }
 \item{collapse_mode}{
	The method to collapse duplicate ids. "mean", "median", "maxSD", "maxIQR", "max" and "min" represent the mean, median, max standard deviation, max interquartile range, maximum and minimum of values for ids in each sample. The default is "maxSD".
  }
}

\value{
	 If there is no format error, the function will return the standardized data matrix and sample annotation data.  Otherwise, the function will output the detailed position of the errors.
}

\note{
	If the users set inputMat as "", the testFileFormat function only test format of sample annotation data. If the users set sampleAnn as "", the testFileFormat function only test format of data matrix. 
}

\author{
	Jing Wang
}

\examples{
	inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
	sampleAnnDir <- system.file("extdata","sampleAnnotation.tsi",package="NetSAM")
	
	formatedData <- testFileFormat(inputMat=inputMatDir,sampleAnn=sampleAnnDir,collapse_mode="maxSD")
}

\keyword{ methods }
