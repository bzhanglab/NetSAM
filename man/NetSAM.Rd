\name{NetSAM}
\alias{NetSAM}

\title{
	Network Seriation and Modularization
}

\description{
	The NetSAM function uses random walk distance-based hierarchical clustering to identify the hierarchical modules of a weighted or unweighted network and then uses the optimal leaf ordering (OLO) method to optimize the one-dimensional ordering of the genes in each module by minimizing the sum of the pair-wise random walk distance of adjacent genes in the ordering. 
}

\usage{
NetSAM(inputNetwork, outputFileName, outputFormat="nsm", edgeType="unweighted", map_to_genesymbol=FALSE, organism="hsapiens", idType="auto",minModule=0.003, stepIte=FALSE, maxStep=4, moduleSigMethod="cutoff", modularityThr=0.2, ZRanNum=10, PerRanNum=100, ranSig=0.05, edgeThr=(-1), nodeThr=(-1), nThreads=3)
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
  \item{outputFormat}{
	The format of the output file. "nsm" format can be used as an input in NetGestalt; "gmt" format can be used to do other network analysis (e.g. as an input in GSEA (Gene Set Enrichment Analysis) to do module enrichment analysis); "multiple" represents the NetSAM function will output three files: ruler file containing gene order information, hmi file containing module information and net file containing network information; and "none" represents the function will not output any file.
  }
  \item{map_to_genesymbol}{
        Because pathway enrichment analysis in NetGestalt is based on gene symbol, setting \code{map_to_genesymbol} as TRUE can transform other ids in the network into gene symbols and thus allow users to do functional analysis based on the identified modules. If the input network is not a biology network or users do not plan to do enrichment analysis in the NetGestalt, users can set \code{map_to_genesymbol} as FALSE. The default is FALSE.
  }
  \item{organism}{
	The organism of the input network. Currently, the package supports the following nine organisms: hsapiens, mmusculus, rnorvegicus, drerio, celegans, scerevisiae, cfamiliaris, dmelanogaster and athaliana. The default is "hsapiens". 
  }
   \item{idType}{
  	The id type of the ids in the input network. MatSAM will use BiomaRt package to transform the input ids to gene symbols based on \code{idType}. The users can also set \code{idType} as "auto" that means MatSAM will automatically search the id type of the input data. However, this may take 10 minutes based on the users' internet speed. The default is "auto".
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
	The significance level for determining whether a network has non-random internal modular organization for the "zscore" or "permutation" methods.
  }
  \item{edgeThr}{
	If the network is too big, it will take a long time to identify the modules. Thus, the function provides the option to set the threshold of number of edges and nodes as \code{edgeThr} and \code{nodeThr}. If the size of network is over the threshold, the function will stop and the users should change the parameters and re-run the function. We suggest to set the threshold for node as 12,000 and the threshold for edge as 300,000. The default is -1, which means there is no limitation for the input network.
  }
  \item{nodeThr}{
  	see \code{edgeThr}.
  }
  \item{nThreads}{
	NetSAM function supports parallel computing based on multiple cores. The default is 3.
  }
}

\value{
	If output format is "nsm", the function will output not only a "nsm" file but also a list object containing module information, gene order information and network information. If output format is "gmt", the function will output the "gmt" file and a matrix object containing the module and annotation information.
}

\note{
	Because the seriation step requires pair-wise distance between all nodes, NetSAM is memory consuming. We recommend to use the 64 bit version of R to run the NetSAM. For networks with less than 10,000 nodes, we recommend to use a computer with 8GB memory. For networks with more than 10,000 nodes, a computer with at least 16GB memory is recommended.
}

\author{
	Jing Wang
}

\examples{
	inputNetworkDir <- system.file("extdata","exampleNetwork.net",package="NetSAM")
	outputFileName <- paste(getwd(),"/NetSAM",sep="")
	result <- NetSAM(inputNetwork=inputNetworkDir, outputFileName=outputFileName, outputFormat="nsm", edgeType="unweighted", map_to_genesymbol=FALSE, organism="hsapiens", idType="auto",minModule=0.003, stepIte=FALSE, maxStep=4, moduleSigMethod="cutoff", modularityThr=0.2, ZRanNum=10, PerRanNum=100, ranSig=0.05, edgeThr=(-1), nodeThr=(-1), nThreads=3)
}

\keyword{ methods }
