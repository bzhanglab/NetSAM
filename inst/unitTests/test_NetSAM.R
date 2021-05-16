test_NetSAM <- function(){
	inputNetworkDir <- system.file("extdata","exampleNetwork.net",package="NetSAM")
	outputFileName <- paste(getwd(),"/NetSAM",sep="")
	result <- NetSAM(inputNetwork=inputNetworkDir, outputFileName=outputFileName)
	checkTrue(!is.na(result[1]))
	checkTrue(!is.na(result[2]))
	checkTrue(!is.na(result[3]))
}
