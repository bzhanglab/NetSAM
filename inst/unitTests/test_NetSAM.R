test_NetSAM <- function(){
	data(inputNetwork)
	outputFileName <- paste(getwd(),"/NetSAM",sep="")
	result <- NetSAM(inputNetwork, outputFileName)
	checkTrue(!is.na(result[1]))
	checkTrue(!is.na(result[2]))
	checkTrue(!is.na(result[3]))
}
