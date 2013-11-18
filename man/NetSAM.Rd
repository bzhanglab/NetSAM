\name{NetSAM}
\alias{NetSAM}

\title{
	Network Seriation and Modularization
}

\description{
	The NetSAM function uses random walk distance-based hierarchical clustering to identify the hierarchical modules of the network and then uses the optimal leaf ordering (OLO) method to optimize the one-dimensional ordering of the genes in each module by minimizing the sum of the pair-wise random walk distance of adjacent genes in the ordering. 
}

\usage{
NetSAM(inputNetwork, outputFileName, minModule = (-1), maxStep = 4, method = "Modularity Cutoff", ModularityThr = 0.2, ZRandomNum = 10, permuteNum = 100, pThr = 0.05)
}

\arguments{
  \item{inputNetwork}{
	The network under analysis. inputNetwork can be the name of the input network file in the edge-list format (each row represents an edge with two node names separated by a tab or space) or be a data object in R (data object must be graphNEL class or data.frame (or matrix) class with two columns).
  }
  \item{outputFileName}{
	The name of the output file. If no file path is provided, the output file will be saved to the current path.
  }
  \item{minModule}{
	The minimum number of nodes for a module (or minimum module size). If the size of a module identified by the function is less than the specified number, the module will not be further partitioned into sub-modules. The default is -1 which means NetSAM will set \code{minModule} as 5 or 0.3 percent of the number of nodes in the input network, whichever is larger.
  }
  \item{maxStep}{
	Because NetSAM uses random walk distance-based hierarchical clustering to reveal the hierarchical organization of an input network, it requires a specified length of the random walks. To get the optimal length, the function will test a range of lengths ranging from 2 to \code{maxStep}. The default is 4.
  }
  \item{method}{
	To test whether a network under consideration has a non-random internal modular organization, the function provides three options: "Modularity Cutoff", "ZScore" and "Permutation". "Modularity Cutoff" means if the modularity score of the network is above a specified cutoff value, the network will be considered to have internal organization and will be further partitioned. For "ZScore" and "Permutation", the function first uses the edge switching method to generate a given number of random networks with the same number of nodes and an identical degree sequence and calculates the modularity scores for these random networks. Then, "ZScore" method will transform the real modularity score to a z score based on the random modularity scores and then transform the z score to a p value assuming a standard normal distribution. The "Permutation" method will compare the real modularity score with the random ones to calculate a p value. Finally, under a specified significance level, the function determines whether the network can be further partitioned. The default is "Modularity Cutoff".
  }
  \item{ModularityThr}{
	Threshold of modularity score for the  "Modularity Cutoff" method. The default is 0.2
  }
  \item{ZRandomNum}{
	The number of random networks that will be generated for the "ZScore" calculation. The default is 10.
  }
  \item{permuteNum}{
	The number of random networks that will be generated for the "Permutation" p value calculation. The default is 100.
  }
  \item{pThr}{
	The significance level for determining whether a network has non-random internal modular organization for the "ZScore" or "Permutation" methods.
  }
}

\note{
	Because the seriation step requires pair-wise distance between all nodes, NetSAM is memory consuming. We recommend to use the 64 bit version of R to run the NetSAM. For networks with less than 10,000 nodes, we recommend to use a computer with 8GB memory. For networks with more than 10,000 nodes, a computer with at least 16GB memory is recommended.
}

\author{
	Jing Wang
}

\examples{
	cat("The input network can be a file.\n")
	inputNetwork <- system.file("extdata","exampleNetwork.txt",package="NetSAM")
	cat("The input network can be also a data object, such as graphNEL object.\n")
	data(inputNetwork)
	outputFileName <- paste(getwd(),"/NetSAM",sep="")
	result <- NetSAM(inputNetwork, outputFileName, minModule = (-1), maxStep = 4, method = "Modularity Cutoff", ModularityThr = 0.2, ZRandomNum = 10, permuteNum = 100, pThr = 0.05)
}

\keyword{ methods }
