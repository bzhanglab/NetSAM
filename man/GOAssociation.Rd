\name{GOAssociation}
\alias{GOAssociation}

\title{
	Identify the associated GO terms for each module
}

\description{
	 The GOAssociation function will identify the associated GO terms for each module from the NetSAM or MatSAM function.
}

\usage{
	GOAssociation(NetSAMOutput, outputHtmlFile, organism, outputType, fdrmethod, fdrth, topNum)
}

\arguments{
 \item{NetSAMOutput}{
	The list object outputted from the NetSAM or MatSAM function.
  }
 \item{outputHtmlFile}{
	The output directory of the HTML file.
  }
 \item{organism}{
	The organism of the input data matrix that has been used to identify the modules. Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana. The default is "hsapiens". 
  }
 \item{outputType}{
	The function supports two types of output results. 1. "significant" represents all associated GO terms should be significant under a certain FDR threshold; 2. "top" represents the function first sorts all GO terms based on their hypergenometric test p values and then selects top GO terms as the associated terms. The default is "significant".
  }
 \item{fdrmethod}{
	The FDR method for identifying the significantly associated GO terms. The default is "BH".
  }
 \item{fdrth}{
	The FDR threshold.
  }
 \item{topNum}{
	The number of the selected top GO terms.
  }

}

\value{
	The function will output a data.frame object and a HTML file to show the associated GO terms for each module.
}

\author{
	Jing Wang
}

\examples{
	data(NetSAMOutput_Example)
	outputHtmlFile <- paste(getwd(),"/GOAsso_HTML",sep="")
	GOAsso <- GOAssociation(NetSAMOutput=netsam_output, outputHtmlFile=outputHtmlFile, organism="hsapiens", fdrmethod="BH", fdrth=0.05, topNum=5)
}

\seealso{
	\code{\link{MatSAM}}
	\code{\link{NetSAM}}
}

\keyword{ methods }
