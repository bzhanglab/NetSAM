\name{NetAnalyzer}
\alias{NetAnalyzer}

\title{
	Network analyzer
}

\description{
	The NetAnalyzer function can calculate the degree, clustering coefficient, betweeness and closeness centrality for each node and the shortest path distance for each pair of nodes. The NetAnalyzer function can also plot the distributions for these measurements.
}

\usage{
NetAnalyzer(inputNetwork, outputFileName,edgeType="unweighted")
}

\arguments{
  \item{inputNetwork}{
	The network under analysis. \code{inputNetwork} can be the directory of the input network file including the file name with "net" extension. If \code{edgeType} is "unweighted", each row represents an edge with two node names separated by a tab or space. If \code{edgeType} is "weighted", each row represents an edge with two node names and edge weight separated by a tab or space.  \code{inputNetwork} can also be a data object in R (data object must be igraph, graphNEL, matrix or data.frame class). 
  }
  \item{edgeType}{
	The type of the input network: "weighted" or "unweighted". 
  }
  \item{outputFileName}{
	The name of the output file.
  }
}

\value{
	The function will output two "txt" files and five "pdf" files. Two "txt" files contain degree, clustering coefficient, betweeness and closeness centrality for each node and the shortest path distance for each pair of nodes. Five "pdf" files are the distributions of these measurements.
}

\author{
	Jing Wang
}

\examples{
	inputNetworkDir <- system.file("extdata","exampleNetwork.net",package="NetSAM")
	outputFileName <- paste(getwd(),"/NetSAM",sep="")
	NetAnalyzer(inputNetworkDir,outputFileName,"unweighted")
}

\keyword{ methods }
