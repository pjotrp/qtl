#####################################################################
#
# PipelineMolgenis.R
#
# copyright (c) 2009, Danny Arends
# last modified Mrt, 2009
# first written Mrt, 2009
# 
# Part of the R/qtl package
# Contains: PipelineMolgenis
#
######################################################################

######################################################################
#
# PipelineMolgenis:
#
######################################################################

PipelineMolgenis <- function(DBmarkerID=298,DBtraitID=181,saveAS="MQMResults",DBpath="http://celtic.service.rug.nl:8080/molgenis4rsandbox",Fupdate=0,...){
	cat("INFO:Starting Stage1: Getting information from Molgenis\n")
		cross <- CrossFromMolgenis(DBmarkerID,DBtraitID,DBpath)
	cat("INFO:DONE Stage1\n")
	cat("INFO:Starting Stage2: scanning for QTL's\n")	
		result <- scanMQMall(cross,...,plot=FALSE)
	cat("INFO:DONE Stage2\n")	
	cat("INFO:Stage3: Storing calculated QTL's to Molgenis\n")
		ResultsToMolgenis(result, saveAS, DBpath, Fupdate)
	cat("INFO:DONE Stage3\n")
}