\name{featureAssociation}
\alias{featureAssociation}

\title{
	Calculate the associations between modules and sample features
}

\description{
	 The featureAssociation function will calculate the associations between the sample features in the input annotation data and the modules identified by NetSAM or MatSAM function.
}

\usage{
	featureAssociation(inputMat, sampleAnn, NetSAMOutput, outputHtmlFile, CONMethod="spearman", CATMethod="kruskal", BINMethod="ranktest", fdrmethod="BH",pth=0.05,collapse_mode="maxSD")
}

\arguments{
  \item{inputMat}{
	\code{inputMat} should contain a file name with extension "cct" or "cbt" or a matrix or data.frame object in R. The first column and first row of the "cct" or "cbt" file should be the row and column names, respectively and other parts are the numeric values. The detail information of "cct" or "cbt" format can be found in the manual of NetGestalt (www.netgestalt.org). A matrix or data.frame object should have row and column names and only contain numeric or integer values.
  }
  \item{sampleAnn}{
	\code{sampleAnn} should be a directory containing a file name with  "tsi" extension or a data.frame object in R. The detail information of "tsi" format can be found in the manual of NetGestalt (www.netgestalt.org). The first row of the sample annotation data is the feature names. The second row is the feature types. The function supports four types: BIN (binary feature, such as male and female), CAT (category feature, such as stage i, stage ii and stage iii), CON (continuous feature, such as age) and SUR (survival data, such as overall survival). The third row is the feature categories. If there is no category information for the features, the sample information will start from the third row . The first column is the sample names.
  }
 \item{NetSAMOutput}{
	The list object outputted from NetSAM or MatSAM.
  }
 \item{outputHtmlFile}{
	The output directory of the HTML file.
  }
 \item{CONMethod}{
	The method to calculate the associations between modules and continuous features. The function provides two methods: "spearman" and "pearson" and the default is "spearman".
  }
 \item{CATMethod}{
	The method to calculate the associations between modules and category features. The function provides two methods: "anova" and "kruskal" and the default is "kruskal".
  }
 \item{BINMethod}{
	The method to calculate the associations between modules and binary features. The function provides two methods: "test" and "rankest" and the default is "rankest".
  }
\item{fdrmethod}{
	The FDR method for identifying the significantly associated GO terms. The default is "BH".
}
 \item{pth}{
	The threshold of the p values to identify the significant associations.
  }
 \item{collapse_mode}{
	The method to collapse duplicate ids. "mean", "median", "maxSD", "maxIQR", "max" and "min" represent the mean, median, max standard deviation, max interquartile range, maximum and minimum of values for ids in each sample. The default is "maxSD".
  }

}

\value{
	The function will output a data.frame object and a HTML file to show the significant associations.
}

\author{
	Jing Wang
}

\examples{
	inputMatDir <- system.file("extdata","exampleExpressionData.cct",package="NetSAM")
	sampleAnnDir <- system.file("extdata","sampleAnnotation.tsi",package="NetSAM")
	data(NetSAMOutput_Example)
	outputHtmlFile <- paste(getwd(),"/featureAsso_HTML",sep="")
	featureAsso <- featureAssociation(inputMat= inputMatDir, sampleAnn=sampleAnnDir, NetSAMOutput=netsam_output, outputHtmlFile=outputHtmlFile, CONMethod="spearman", CATMethod="kruskal", BINMethod="ranktest", fdrmethod="BH",pth=0.05,collapse_mode="maxSD")
}

\seealso{
	\code{\link{MatSAM}}
	\code{\link{NetSAM}}
}

\keyword{ methods }
